[Mesh]
  type = GeneratedMesh
  dim = 1
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [null]
    type = NullKernel
    variable = u
  []
[]

[OptimizationReporter]
  type = VariationalInferenceOptimization
  objective_name = objective_value
  parameter_names = 'p1 p2'
  num_values = '2 3'
  initial_condition = '1 2; 4 5 6'
  upper_bounds = '100; 200'
  lower_bounds = '-1; -2'
  num_measurements_per_experiment = 1000
  experimental_noise = 11610.00
  num_experiments = 4
  elbo_objective_name = 'ELBO'
  outputs = out
[]

[UserObjects]
  [optReporterTester]
    type = OptimizationReporterTest
    values_to_set_parameters_to = '100 100 200 200 200 1 0 1 1 0 1 0 0 1'
    expected_lower_bounds = '-1 -1 -1 -1 -1 -2 -2 -2 -2 -2 -2 -2 -2 -2'
    expected_upper_bounds = '100 100 100 100 100 200 200 200 200 200 200 200 200 200'
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  [out]
    type = JSON
    execute_system_information_on = none
  []
[]
