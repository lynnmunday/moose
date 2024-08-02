[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = 10
    nx = 10
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[AuxVariables]
  [u]
  []
  [du_dp1]
  []
  [du_dp2]
  []
[]
[AuxKernels]
  [u]
    type = FunctorAux
    variable = u
    functor = u_func
  []
  [du_dp1]
    type = FunctorAux
    variable = du_dp1
    functor = du_dp1_func
  []
  [du_dp2]
    type = FunctorAux
    variable = du_dp2
    functor = du_dp2_func
  []
[]

[Functions]
  #u is a "fake" FEM simulation result
  [u_func]
    type = ParsedOptimizationFunction
    expression = 'p1*sin(x)+p2*cos(x)'
    param_symbol_names = 'p1 p2'
    param_vector_name = vals/vals
  []
  [du_dp1_func]
    type = ParsedFunction
    expression = 'sin(x)'
  []
  [du_dp2_func]
    type = ParsedFunction
    expression = 'cos(x)'
  []
[]

[Reporters]
  #objective values
  [u_m]
    type = OptimizationData
    variable = u
    measurement_file = 'allExpMeasurements.csv'
    file_xcoord = 'x'
    file_ycoord = 'y'
    file_zcoord = 'z'
    file_value = 'value'
  []
  # J = 1/2||u-u_m||
  # u = simulation value
  # u_m measurement values for all experiments
  [objective]
    type = ParsedVectorRealReductionReporter
    name = value
    reporter_name = u_m/misfit_values
    initial_value = 0
    expression = 'reduction_value+(0.5*indexed_value*indexed_value)'
    outputs = json
  []

  #gradient values
  # this is the average value over all experiments at each measurement point
  # need avg_misfit = (u-u_m_avg)^2
  [u_m_avg]
    type = OptimizationData
    variable = u
    measurement_values = '1.54024167 0.97013352 -0.49004537 -1.50395413 -1.13178525 0.27719184'
    measurement_points = '1 0 0
    2 0 0
    3 0 0
    4 0 0
    5 0 0
    6 0 0'
  []

  [L2_dudp1]
    type = ParsedVectorReporter
    name = vec
    reporter_names = 'u_m_avg/misfit_values sample_grad_u/du_dp1'
    reporter_symbols = 'avg_misfit dudp1'
    expression = '-avg_misfit*avg_misfit*dudp1'
  []
  [L2_dudp2]
    type = ParsedVectorReporter
    name = vec
    reporter_names = 'u_m_avg/misfit_values sample_grad_u/du_dp2'
    reporter_symbols = 'avg_misfit dudp2'
    expression = '-avg_misfit*avg_misfit*dudp2'
  []
  [L2_dudp1_sum]
    type = ParsedVectorRealReductionReporter
    name = vec
    reporter_name = L2_dudp1/vec
    initial_value = 0
    expression = 'reduction_value+indexed_value'
  []
  [L2_dudp2_sum]
    type = ParsedVectorRealReductionReporter
    name = vec
    reporter_name = L2_dudp2/vec
    initial_value = 0
    expression = 'reduction_value+indexed_value'
  []
[]

[VectorPostprocessors]
  [sample_grad_u]
    type = PointValueSampler
    variable = 'du_dp1 du_dp2'
    points = '1 0 0
      2 0 0
      3 0 0
      4 0 0
      5 0 0
      6 0 0'
    sort_by = x
    outputs = none
  []
  [grad_f]
    type = VectorOfPostprocessors
    reporters = 'L2_dudp1_sum/vec L2_dudp2_sum/vec'
  []
[]

[Reporters]
  [vals]
    type = ConstantReporter
    real_vector_names = 'vals'
    real_vector_values = '1 0.5'
  []
[]

[Outputs]
  csv = true
  [json]
    type = JSON
    execute_system_information_on = NONE
  []
[]
