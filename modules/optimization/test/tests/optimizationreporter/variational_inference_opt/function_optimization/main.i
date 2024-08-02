[Optimization]
[]

[OptimizationReporter]
  type = VariationalInferenceOptimization
  objective_name = 'obj_value'
  elbo_objective_name = 'ELBO'
  parameter_names = 'vals'
  num_values = '2'
  initial_condition = 1
  num_measurements_per_experiment = 1000
  experimental_noise = 5.00
  num_experiments = 4
  outputs = out
[]

[Problem]
  solve = false
[]
[Executioner]
  type = Optimize
  # tao_solver = taobqnls
  # petsc_options_iname = '-tao_gatol -tao_ls_type'
  # petsc_options_value = '1e-8 unit'

  # tao_solver = taobncg
  # petsc_options_iname = '-tao_max_it -tao_gatol -tao_ls_type'
  # petsc_options_value = '1 1e-8 unit'

  ## THESE OPTIONS ARE FOR TESTING THE ADJOINT GRADIENT
  tao_solver = taobncg
  petsc_options_iname = '-tao_max_it -tao_fd_test -tao_test_gradient -tao_fd_gradient -tao_fd_delta -tao_gatol -tao_ls_type'
  petsc_options_value = '1 true true false 1e-10 0.1 unit'
  petsc_options = '-tao_test_gradient_view'

  verbose = true
  output_optimization_iterations = true
[]

[MultiApps]
  [forward]
    type = FullSolveMultiApp
    input_files = forward.i
    execute_on = FORWARD
  []
[]

[Transfers]
  [toForward]
    type = MultiAppReporterTransfer
    to_multi_app = forward
    from_reporters = 'OptimizationReporter/vals'
    to_reporters = 'vals/vals'
  []

  [fromForward]
    type = MultiAppReporterTransfer
    from_multi_app = forward
    from_reporters = 'objective/value
                      grad_f/grad_f'
    to_reporters = 'OptimizationReporter/obj_value
                    OptimizationReporter/grad_vals'
  []
[]

[Outputs]
  [json]
    type = JSON
    execute_system_information_on = none
  []
  [json_forward]
    type = JSON
    execute_on = 'FORWARD'
    execute_system_information_on = none
  []
[]
