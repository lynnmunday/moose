#include "TabulatedPolynomialModel3.h"
#include <fstream>

namespace neml2
{
register_NEML2_object(TabulatedPolynomialModel3);

OptionSet
TabulatedPolynomialModel3::expected_options()
{
  auto options = Model::expected_options();
  // Model inputs
  options.set<VariableName>("von_mises_stress") = VariableName("state", "s");
  options.set<VariableName>("temperature") = VariableName("forces", "T");
  options.set<VariableName>("equivalent_plastic_strain") = VariableName("state", "ep");
  options.set<VariableName>("internal_state_1") = VariableName("state", "s1");
  options.set<VariableName>("internal_state_2") = VariableName("state", "s2");
  options.set<VariableName>("internal_state_3") = VariableName("state", "s3");
  // Model outputs
  options.set<VariableName>("equivalent_plastic_strain_rate") = VariableName("state", "ep_rate");
  options.set<VariableName>("internal_state_1_rate") = VariableName("state", "s1_rate");
  options.set<VariableName>("internal_state_2_rate") = VariableName("state", "s2_rate");
  // JSON
  options.set<std::string>("model_file_name"); // = "../../../data/laromance/test/3tile.json";
  // Use AD
  options.set<bool>("_use_AD_first_derivative") = true;
  options.set<bool>("_use_AD_second_derivative") = true;
  return options;
}

TabulatedPolynomialModel3::TabulatedPolynomialModel3(const OptionSet & options)
  : Model(options),
    _s(declare_input_variable<Scalar>("von_mises_stress")),
    _T(declare_input_variable<Scalar>("temperature")),
    _ep(declare_input_variable<Scalar>("equivalent_plastic_strain")),
    _s1(declare_input_variable<Scalar>("internal_state_1")),
    _s2(declare_input_variable<Scalar>("internal_state_2")),
    _s3(declare_input_variable<Scalar>("internal_state_3")),
    _ep_dot(declare_output_variable<Scalar>("equivalent_plastic_strain_rate")),
    _s1_dot(declare_output_variable<Scalar>("internal_state_1_rate")),
    _s2_dot(declare_output_variable<Scalar>("internal_state_2_rate"))
{
  std::string filename = options.get<std::string>("model_file_name");
  std::ifstream model_file(filename.c_str());
  model_file >> _json;

  // storing grid points for indexing.
  // these should be stored differently so that they are all read in at once.  The order of this can
  // get messed up easily
  _stress_grid = json_vector_to_torch("in_stress");
  _temperature_grid = json_vector_to_torch("in_temperature");
  _plastic_strain_grid = json_vector_to_torch("in_plastic_strain");
  _cell_grid = json_vector_to_torch("in_cell");
  _wall_grid = json_vector_to_torch("in_wall");
  _env_grid = json_vector_to_torch("in_env");

  // Storing values for interpolation
  _grid_values_ep = json_6Dvector_to_torch("out_ep");
  _grid_values_cell = json_6Dvector_to_torch("out_cell");
  _grid_values_wall = json_6Dvector_to_torch("out_wall");
}

void
TabulatedPolynomialModel3::set_value(bool out, bool dout_din, bool d2out_din2)
{
  neml_assert_dbg(!dout_din || !d2out_din2,
                  "Only AD derivatives are currently supported for this model");

  if (out)
  {
    //"in_stress": [0,...,300]
    //"in_temperature": [600,...,1100]
    // "in_plastic_strain": [0,...,0.01]
    // "in_cell": [1485051.0319252764, 8500000000000.0]
    // "in_wall": [5033298671083.821, 12000000000000.0]
    // "in_env": [9.99999999999973e-10, 1e-06]
    // I think this or the function findLeftIndexAndFraction needs to have a type for scalar so that
    // it is still a klong and not double
    std::vector<std::pair<Scalar, Scalar>> left_index_weight;
    left_index_weight.push_back(findLeftIndexAndFraction(_stress_grid, _s));
    left_index_weight.push_back(findLeftIndexAndFraction(_temperature_grid, _T));
    left_index_weight.push_back(findLeftIndexAndFraction(_plastic_strain_grid, _ep));
    left_index_weight.push_back(findLeftIndexAndFraction(_cell_grid, _s1));
    left_index_weight.push_back(findLeftIndexAndFraction(_wall_grid, _s2));
    left_index_weight.push_back(findLeftIndexAndFraction(_env_grid, _s3));

    _ep_dot = compute_interpolation(left_index_weight, _grid_values_ep);
    _s1_dot = compute_interpolation(left_index_weight, _grid_values_cell);
    _s2_dot = compute_interpolation(left_index_weight, _grid_values_wall);
    // these need transformed and turned into rates or whatever should be returned
  }
}

std::pair<Scalar, Scalar>
TabulatedPolynomialModel3::findLeftIndexAndFraction(const Scalar grid, const Scalar interp_points)
{
  // idx is for the left grid point.
  // searchsorted returns the right idx so -1 makes it the left
  Scalar left_idx =
      Scalar(torch::searchsorted(grid, interp_points)) - torch::tensor(1, torch::kLong);

  Scalar left_coord = grid.index({left_idx.to(torch::kLong)});
  Scalar right_coord = grid.index({left_idx.to(torch::kLong) + torch::tensor(1, torch::kLong)});
  // fixme lynn, I should be able to do this above but I can't get it to work in the
  // left_coord fraction is for the left fraction.  For the above example, fraction
  // would be 0.8--> 0.8*20+0.2*30=22
  Scalar left_fraction = (right_coord - interp_points) / (right_coord - left_coord);
  return std::pair<Scalar, Scalar>(left_idx, torch::stack({left_fraction, 1 - left_fraction}, -1));
}

Scalar
TabulatedPolynomialModel3::compute_interpolation(
    const std::vector<std::pair<Scalar, Scalar>> index_and_fraction, const Scalar grid_values)
{
  // need to get batch size from variable
  Scalar result = Scalar::zeros_like(_T);
  // Manually unrolled loops for 6D
  for (auto i = 0; i < 2; ++i)
    for (auto j = 0; j < 2; ++j)
      for (auto k = 0; k < 2; ++k)
        for (auto l = 0; l < 2; ++l)
          for (auto m = 0; m < 2; ++m)
            for (auto n = 0; n < 2; ++n)
            {
              // I shouldn't have to cast these to klong,  above left_fraction is klong but loses
              // the type when it is returned from the function
              auto vertex_value =
                  grid_values.index({(index_and_fraction[0].first + i).to(torch::kLong),
                                     (index_and_fraction[1].first + j).to(torch::kLong),
                                     (index_and_fraction[2].first + k).to(torch::kLong),
                                     (index_and_fraction[3].first + l).to(torch::kLong),
                                     (index_and_fraction[4].first + m).to(torch::kLong),
                                     (index_and_fraction[5].first + n).to(torch::kLong)});
              auto weight = index_and_fraction[0].second.select(-1, i) *
                            index_and_fraction[1].second.select(-1, j) *
                            index_and_fraction[2].second.select(-1, k) *
                            index_and_fraction[3].second.select(-1, l) *
                            index_and_fraction[4].second.select(-1, m) *
                            index_and_fraction[5].second.select(-1, n);
              result += vertex_value * weight;
            }
  return result;
}

Scalar
TabulatedPolynomialModel3::json_vector_to_torch(std::string key)
{
  std::cout << key << std::endl;
  neml_assert_dbg(_json.contains(key), "The key '", key, "' is missing from the JSON data file.");
  std::vector<Real> in_data = _json[key].template get<std::vector<Real>>();
  const long long data_dim = in_data.size();
  return Scalar(torch::from_blob(in_data.data(), {data_dim}).clone());
}

Scalar
TabulatedPolynomialModel3::json_6Dvector_to_torch(std::string key)
{
  std::cout << key << std::endl;
  neml_assert_dbg(_json.contains(key), "The key '", key, "' is missing from the JSON data file.");

  std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Real>>>>>> out_data =
      _json[key]
          .template get<
              std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Real>>>>>>>();

  std::vector<Real> linearize_values;
  for (auto && level0 : out_data)
    for (auto && level1 : level0)
      for (auto && level2 : level1)
        for (auto && level3 : level2)
          for (auto && level4 : level3)
            for (auto && value : level4)
              linearize_values.push_back(value);

  long long sz_l0 = out_data.size();
  long long sz_l1 = out_data[0].size();
  long long sz_l2 = out_data[0][0].size();
  long long sz_l3 = out_data[0][0][0].size();
  long long sz_l4 = out_data[0][0][0][0].size();
  long long sz_l5 = out_data[0][0][0][0][0].size();
  return Scalar(
      torch::from_blob(linearize_values.data(), {sz_l0, sz_l1, sz_l2, sz_l3, sz_l4, sz_l5})
          .clone());
}
} // namespace neml2
