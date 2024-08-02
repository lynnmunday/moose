//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseError.h"
#include "MooseTypes.h"
#include "VariationalInferenceOptimization.h"
#include "libmesh/id_types.h"
#include <numeric>
#include <random>
registerMooseObject("OptimizationApp", VariationalInferenceOptimization);

InputParameters
VariationalInferenceOptimization::validParams()
{
  InputParameters params = GeneralOptimization::validParams();

  params.addRequiredParam<ReporterValueName>(
      "objective_name", "Preferred name of reporter value defining the objective.");
  params.addParam<std::vector<dof_id_type>>(
      "num_values",
      "Number of parameter values associated with each parameter group in 'parameter_names'.");
  params.addParam<ReporterValueName>("num_values_name",
                                     "Reporter that holds the number of parameter values "
                                     "associated with each parameter group in 'parameter_names'.");
  params.addParam<std::vector<std::vector<Real>>>(
      "initial_condition",
      "Initial conditions for each parameter. A vector is given for each parameter group.  A "
      "single value can be given for each group and all parameters in that group will be set to "
      "that value.  The default value is 0.");
  params.addParam<std::vector<std::vector<Real>>>(
      "lower_bounds",
      "Lower bound for each parameter.  A vector is given for each parameter group.  A single "
      "value can be given for each group and all parameters in that group will be set to that "
      "value");
  params.addParam<std::vector<std::vector<Real>>>(
      "upper_bounds",
      "Upper bound for each parameter.  A vector is given for each parameter group.  A single "
      "value can be given for each group and all parameters in that group will be set to that "
      "value");
  params.addParam<unsigned int>("seed", 0, "Seed for random number generator.");

  params.addParam<Real>("experimental_noise", 1, " Mean sample variance of experiments");
  params.addParam<unsigned int>("num_experiments", 1, "Number of experiments");
  params.addRequiredParam<unsigned int>("num_measurements_per_experiment",
                                        "Number of measurement points in each experiment");
  params.addRequiredParam<ReporterValueName>(
      "elbo_objective_name", "Preferred name of reporter value defining the ELBO objective.");
  params.addClassDescription("Reporter that provides TAO with the objective, gradient, and "
                             "constraint data, which are supplied by the reporters and "
                             "postprocessors from the forward and adjoint subapps.");

  return params;
}

VariationalInferenceOptimization::VariationalInferenceOptimization(
    const InputParameters & parameters)
  : GeneralOptimization(parameters),
    _elbo_objective_val(declareValueByName<Real>(getParam<ReporterValueName>("elbo_objective_name"),
                                                 REPORTER_MODE_REPLICATED)),
    _seed(getParam<unsigned int>("seed")),
    _experimental_noise(getParam<Real>("experimental_noise")),
    _num_experiments(getParam<unsigned int>("num_experiments")),
    _num_measurements_per_experiment(getParam<unsigned int>("num_measurements_per_experiment"))
{
  for (const auto & i : make_range(_nparams))
  {
    _parameters.push_back(&declareValueByName<std::vector<Real>>(
        "covariance_" + _parameter_names[i], REPORTER_MODE_REPLICATED));
    _sampled_parameters.push_back(&declareValueByName<std::vector<Real>>(
        "sampled_" + _parameter_names[i], REPORTER_MODE_REPLICATED));
    _gradients.push_back(&declareValueByName<std::vector<Real>>(
        "grad_covariance_" + _parameter_names[i], REPORTER_MODE_REPLICATED));
  }
}

Real
VariationalInferenceOptimization::computeObjective()
{

  Real regularization = 0;
  if (_tikhonov_coeff > 0.0)
  {
    Real param_norm_sqr = 0;
    for (const auto & data : _parameters)
      for (const auto & param_val : *data)
        param_norm_sqr += param_val * param_val;
    // We multiply by 0.5 to maintain  backwards compatibility.
    regularization = 0.5 * _tikhonov_coeff * param_norm_sqr;
  }

  // variational inference objective terms
  std::size_t n_params = _parameter_names.size();
  // covariance objective terms
  Real obj_entropy = 0;
  for (std::size_t i = 0; i < n_params; ++i)
  {
    DenseMatrix<Real> L(_group_covariance[i]);
    obj_entropy += std::log(L.det()) + 0.5 * _nvalues[i] * (1 + std::log(2 * M_PI));
  }

  // compute objective terms that only rely on sampled parameters
  Real dot_product = 0;
  for (std::size_t ii = 0; ii < _initial_mean_parameters.size(); ++ii)
    for (std::size_t jj = 0; jj < _initial_mean_parameters[ii].size(); ++jj)
      dot_product += (_sampled_parameters[ii]->at(jj) - _initial_mean_parameters[ii][jj]) *
                     (_sampled_parameters[ii]->at(jj) - _initial_mean_parameters[ii][jj]);

  dot_product = -1.0 * std::log(2 * M_PI) - .5 * dot_product;
  Real likelihood_free_terms = dot_product - obj_entropy;

  // Constant LogLiklihood for objective function
  Real const_term_loglikelihood = -0.5 * (Real)_num_experiments *
                                  (Real)_num_measurements_per_experiment *
                                  std::log(2.0 * M_PI * _experimental_noise);
  // combine and compute ELBO
  _elbo_objective_val = likelihood_free_terms + const_term_loglikelihood -
                        (_objective_val / _experimental_noise + regularization);

  return _elbo_objective_val;
}

void
VariationalInferenceOptimization::computeGradient(libMesh::PetscVector<Number> & gradient) const
{
  // fixme LYNN  _nvalues contains [mean vals,covariance_vals]
  //-- make sure grad is being filled in the correct order
  for (const auto & param_group_id : make_range(_nparams))
  {
    if (_gradients[param_group_id]->size() != _nvalues[param_group_id])
      mooseError("The gradient for parameter ",
                 _parameter_names[param_group_id],
                 " has changed, expected ",
                 _nvalues[param_group_id],
                 " versus ",
                 _gradients[param_group_id]->size(),
                 ".");

    if (_tikhonov_coeff > 0.0)
    {
      auto params = _parameters[param_group_id];
      auto grads = _gradients[param_group_id];
      for (const auto & param_id : make_range(_nvalues[param_group_id]))
        (*grads)[param_id] += (*params)[param_id] * _tikhonov_coeff;
    }
  }

  // modifying for variational inference

  // ******* grad wrt mu terms **********
  // grad mu of log-likelihood
  for (const auto & param_group_id : make_range(_nparams))
  {
    auto grad_all = _gradients[param_group_id];
    for (const auto & param_id : make_range(_nvalues[param_group_id]))
      (*grad_all)[param_id] += (*grad_all)[param_id] / _experimental_noise;
  }

  // grad mu of prior
  std::vector<std::vector<Real>> grad_mu_prior;
  for (const auto & param_group_id : make_range(_nparams))
  {
    std::vector<Real> grad_group;
    for (std::size_t i = 0; i < _initial_mean_parameters[param_group_id].size(); ++i)
    {
      grad_group.push_back(_initial_mean_parameters[param_group_id][i] -
                           _sampled_parameters[param_group_id]->at(i));
    }
    grad_mu_prior.push_back(grad_group);
  }
  // add grad mu of prior into mean terms in full gradient
  for (const auto & param_group_id : make_range(_nparams))
  {
    auto grad_all = _gradients[param_group_id];
    auto grad_mu = grad_mu_prior[param_group_id];
    for (const auto & param_id : make_range(_nvalues[param_group_id]))
      (*grad_all)[param_id] += grad_mu[param_id];
  }

  // ******* grad wrt L terms **********
  // grad L of log Likelihood
  // this needs grads from forward problem
  std::vector<std::vector<Real>> grad_LL_L;
  for (const auto & param_group_id : make_range(_nparams))
  {
    std::vector<Real> LL_cov_vec;
    DenseMatrix<Real> LL_cov;
    DenseVector<Real> eps(_random_draws[param_group_id]);
    DenseVector<Real> dR_dp(*_gradients[param_group_id]);
    LL_cov.outer_product(dR_dp, eps);
    for (std::size_t j = 0; j < _nvalues[param_group_id]; ++j)
      for (std::size_t k = 0; k < _nvalues[param_group_id]; ++k)
        if (k < j)
          LL_cov_vec.push_back(LL_cov(j, k) / _experimental_noise);
    grad_LL_L.push_back(LL_cov_vec);
  }
  // grad L of prior
  std::vector<std::vector<Real>> grad_prior_L;
  for (const auto & param_group_id : make_range(_nparams))
  {
    DenseVector<Real> mup(_initial_mean_parameters[param_group_id]);
    DenseVector<Real> eps(_random_draws[param_group_id]);
    DenseMatrix<Real> L(_group_covariance[param_group_id]);
    DenseVector<Real> mu(*_parameters[param_group_id]);
    DenseMatrix<Real> mup_eps_T;
    DenseMatrix<Real> L_eps_eps_T;
    DenseMatrix<Real> mu_eps_T;
    mup_eps_T.outer_product(mup, eps);
    L_eps_eps_T.outer_product(eps, eps);
    L_eps_eps_T.left_multiply(L);
    mu_eps_T.outer_product(mu, eps);
    // Set elememts above the diagonal to 0
    // creates a lower-triangular matrix
    for (std::size_t i = 0; i < _nvalues[param_group_id]; ++i)
      for (std::size_t j = 0; j < _nvalues[param_group_id]; ++j)
        if (j > i)
        {
          mup_eps_T(i, j) = 0;
          L_eps_eps_T(i, j) = 0;
          mu_eps_T(i, j) = 0;
        }
    DenseMatrix<Real> dP_dL(mup_eps_T);
    dP_dL -= L_eps_eps_T;
    dP_dL -= mu_eps_T;
    std::vector<Real> cov_storage(_nvalues[param_group_id + _nparams]); // st
    int counter = 0;
    for (std::size_t row_index = 0; row_index < _nvalues[param_group_id]; row_index++)
      for (std::size_t col_index = 0; col_index <= row_index; col_index++)
      {
        cov_storage[counter] = dP_dL(row_index, col_index);
        counter++; // So that you can index into cov_vec[counter]
      }
    grad_prior_L.push_back(cov_storage);
  }

  //  grad L of entropy
  std::vector<std::vector<Real>> grad_entropy_L;
  for (const auto & param_group_id : make_range(_nparams))
  {
    DenseMatrix<Real> L(_group_covariance[param_group_id]);
    std::vector<Real> inverse_vals;
    for (std::size_t i = 0; i < _nvalues[param_group_id]; ++i)
      for (std::size_t j = 0; j <= i; ++j)
      {
        if (i == j)
          inverse_vals.push_back(1.0 / L(i, i));
        else
          inverse_vals.push_back(0.0);
      }
    grad_entropy_L.push_back(inverse_vals);
  }

  // add grad mu terms into covariance terms of full gradient
  for (std::size_t i = _nparams; i < 2 * _nparams; ++i)
  {
    auto grad_all = _gradients[i];
    auto grad_1 = grad_LL_L[i];
    auto grad_2 = grad_prior_L[i];
    auto grad_3 = grad_entropy_L[i];
    mooseAssert(grad_1.size() == grad_2.size(), "gradients must be same size.");
    mooseAssert(grad_1.size() == grad_3.size(), "gradients must be same size.");
    for (const auto & param_id : make_range(_nvalues[i]))
      (*grad_all)[param_id] += (grad_1[param_id] + grad_2[param_id] + grad_3[param_id]);
  }

  // Now copy gradients into petsc vector for TAO
  OptUtils::copyReporterIntoPetscVector(_gradients, gradient);
}

void
VariationalInferenceOptimization::updateParameters(const libMesh::PetscVector<Number> & x)
{
  OptUtils::copyPetscVectorIntoReporter(x, _parameters);
  // increment seed to get new set of parameters
  _random_draws.clear();

  unsigned int new_seed = _seed + _its;
  std::mt19937 gen(new_seed); // Mersenne Twister generator seeded with rd()
  std::normal_distribution<Real> d(0, 1);
  std::size_t n_params = _parameter_names.size();
  for (std::size_t i = 0; i < n_params; ++i)
  {
    std::vector<Real> random_set_draws;
    for (std::size_t j = 0; j < _nvalues[i]; ++j)
      random_set_draws.push_back(d(gen));
    _random_draws.push_back(random_set_draws);
  }

  _group_covariance.clear();
  for (std::size_t i = 0; i < n_params; ++i)
  {
    std::vector<Real> cov_vec = *_parameters[i + n_params];
    DenseMatrix<Real> L(_nvalues[i], _nvalues[i]);
    L.zero();
    int counter = 0;
    for (std::size_t row_index = 0; row_index < _nvalues[i]; row_index++)
      for (std::size_t col_index = 0; col_index <= row_index; col_index++)
      {
        L(row_index, col_index) = cov_vec[counter];
        counter++; // So that you can index into cov_vec[counter]
      }
    _group_covariance.push_back(L);
  }

  for (std::size_t i = 0; i < n_params; ++i)
  {
    DenseVector<Real> mu(_nvalues[i]);
    mu = *_parameters[i];
    std::vector<Real> cov_vec = *_parameters[i + n_params];
    DenseMatrix<Real> L(_group_covariance[i]);
    DenseVector<Real> random_set_draws(_nvalues[i]);
    random_set_draws = _random_draws[i];

    // make this member data later!! Computing grad of entropy wrto L
    std::vector<Real> inverse_vals;
    for (std::size_t i = 0; i < _nvalues[i]; i++)
    {
      inverse_vals.push_back(1.0 / L(i, i)); // FIUOPWHJAOIFHWNAO{FGI:JEWSD{OGSIO }} Change later
    }

    DenseVector<Real> LE(_nvalues[i]);
    L.vector_mult(LE, random_set_draws);
    LE.add(1.0, mu);

    _sampled_parameters[i]->assign(LE.get_values().begin(), LE.get_values().end());
  }
}

void
VariationalInferenceOptimization::setICsandBounds()
{
  // Check that one and only one set of parameters are given.
  // Set here because some derived reporters use a different method of
  // determining numbers of dofs
  if (!(isParamValid("num_values_name") ^ isParamValid("num_values")))
    paramError("Need to supply one and only one of num_values_name or num_values.");

  if (_num_values_reporter)
    _nvalues = *_num_values_reporter;
  else
    _nvalues = getParam<std::vector<dof_id_type>>("num_values");

  // size checks
  if (_parameter_names.size() != _nvalues.size())
    paramError(
        "num_parameters",
        "There should be a number in \'num_parameters\' for each name in \'parameter_names\'.");

  std::vector<Real> initial_conditions_mean_vec(fillParamsVector("initial_condition", 0));
  _lower_bounds = fillParamsVector("lower_bounds", std::numeric_limits<Real>::lowest());
  _upper_bounds = fillParamsVector("upper_bounds", std::numeric_limits<Real>::max());

  // initialize parameters with initial conditions
  std::size_t stride = 0;
  for (std::size_t i = 0; i < _nvalues.size(); ++i)
  {
    _initial_mean_parameters.push_back(
        {initial_conditions_mean_vec.begin() + stride,
         initial_conditions_mean_vec.begin() + stride + _nvalues[i]});
    _sampled_parameters[i]->assign(initial_conditions_mean_vec.begin() + stride,
                                   initial_conditions_mean_vec.begin() + stride + _nvalues[i]);
    stride += _nvalues[i];
  }

  // Now we need to add in the sizes of the covariance matrix for each parameter
  std::size_t num_mean_values = _nvalues.size();
  for (std::size_t i = 0; i < num_mean_values; ++i)
  {
    dof_id_type n = _nvalues[i];
    dof_id_type entries_in_L_covariance_matrix = (n * n + n) / 2;
    _nvalues.push_back(entries_in_L_covariance_matrix);
  }

  std::vector<Real> initial_conditions_covariance = fillInitialCovarianceParamsVector();
  // merge all initial conditions into a single vector
  std::vector<Real> initial_conditions(initial_conditions_mean_vec);
  initial_conditions.insert(initial_conditions.end(),
                            initial_conditions_covariance.begin(),
                            initial_conditions_covariance.end());

  fillBoundsCovarianceParamsVector(_upper_bounds, std::numeric_limits<Real>::max());
  fillBoundsCovarianceParamsVector(_lower_bounds, std::numeric_limits<Real>::lowest());

  // Now update the total size of the optimization system
  _ndof = std::accumulate(_nvalues.begin(), _nvalues.end(), 0);

  stride = 0;
  for (std::size_t i = 0; i < _nvalues.size(); ++i)
  {
    _gradients[i]->resize(_nvalues[i]);
    _parameters[i]->assign(initial_conditions.begin() + stride,
                           initial_conditions.begin() + stride + _nvalues[i]);
    stride += _nvalues[i];
  }
}

void
VariationalInferenceOptimization::fillBoundsCovarianceParamsVector(std::vector<Real> & data_vec,
                                                                   Real value)
{
  // fill with default values
  std::vector<std::vector<Real>> covariance_data;
  for (std::size_t i = _parameter_names.size(); i < _nvalues.size(); ++i)
  {
    dof_id_type n_cov_params = _nvalues[i];
    std::vector<Real> cov_params(n_cov_params, value);
    covariance_data.push_back(cov_params);
  }

  // flatten into single vector
  for (const auto & vec : covariance_data)
    data_vec.insert(data_vec.end(), vec.begin(), vec.end());
}

std::vector<Real>
VariationalInferenceOptimization::fillInitialCovarianceParamsVector()
{
  // fill with default values
  std::vector<Real> covariance_data;
  for (std::size_t i = _parameter_names.size(); i < _nvalues.size(); ++i)
  {
    dof_id_type n_cov_params = _nvalues[i];
    dof_id_type n_cov_matrix_dim = floor(std::sqrt(2 * n_cov_params));
    std::vector<Real> cov_params(n_cov_params, 0);
    std::size_t id_diag = 0;
    for (std::size_t i_row = 0; i_row < n_cov_matrix_dim; ++i_row)
    {
      cov_params[id_diag] = 1;
      id_diag += i_row + 2;
    }
    covariance_data.insert(covariance_data.end(), cov_params.begin(), cov_params.end());
  }
  return covariance_data;
}
