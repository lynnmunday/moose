[Tensors]
  [A0]
    type = Tensor
    values = "1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6
              1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6"
    base_shape = '(2,3,3)'
  []
  [A1]
    type = Tensor
    values = "1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6
              1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6
              1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6
              1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6
              1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6
              1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6"
    base_shape = '(2,3,3,4)'
  []
  [A2]
    type = Tensor
    values = "1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6
              1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6
              1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6
              1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6
              1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6
              1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6 1e-6 2e-6 3e-6 4e-6"
    base_shape = '(2,3,3,4)'
  []
  [s_lb]
    type = Tensor
    values = '0 50'
    base_shape = (2)
  []
  [s_ub]
    type = Tensor
    values = '50 100'
    base_shape = (2)
  []
  [T_lb]
    type = Tensor
    values = '0 300 600'
    base_shape = (3)
  []
  [T_ub]
    type = Tensor
    values = '300 600 1000'
    base_shape = (3)
  []
[]

[Models]
  [rom]
    type = TabulatedPolynomialModel2
    model_file_name = 'models/polyrom.json'
    von_mises_stress = 'state/s'
    temperature = 'forces/T'
    internal_state_1 = 'state/s1'
    internal_state_2 = 'state/s2'
    equivalent_plastic_strain_rate = 'state/ep_rate'
    internal_state_1_rate = 'state/s1_rate'
    internal_state_2_rate = 'state/s2_rate'
    A2 = 'A2'
    stress_tile_lower_bounds = 's_lb'
    stress_tile_upper_bounds = 's_ub'
    temperature_tile_lower_bounds = 'T_lb'
    temperature_tile_upper_bounds = 'T_ub'
    # use_AD_first_derivative = true
    # use_AD_second_derivative = true
  []
[]
