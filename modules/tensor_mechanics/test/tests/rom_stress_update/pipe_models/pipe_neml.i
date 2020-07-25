#This runs with blackbear

[GlobalParams]
  order = SECOND
  family = LAGRANGE
  displacements = 'disp_x disp_y'
[]

[Problem]
  coord_type = RZ
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Mesh]
  [pipe]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 1.884
    xmax = 2.084
    ymin = 0
    ymax = 0.04
    nx = 5
    ny = 1
    elem_type = quad9
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        strain = FINITE
        add_variables = true
        generate_output = 'max_principal_stress max_principal_strain vonmises_stress
                          strain_xx strain_yy strain_zz
                          stress_xx stress_yy stress_zz'
      []
    []
  []
[]

[AuxVariables]
  [saved_x]
  []
  [saved_y]
  []
  [sint]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  [inner_pressure]
   type = PiecewiseLinear
   x = '0 5000 2e8'  #time
#   y = '0.0 1.0e6 1.0e6'  #pressure in Pa
   y = '0.0 1.0 1.0'  #pressure in MPa
   scale_factor=5
  []
  [timefunc]
    type=ParsedFunction
    value = t
  []
[]

[AuxKernels]
  [./sint]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = sint
    scalar_type = stressIntensity
    execute_on = timestep_end
  [../]
[]

[BCs]
#  [fixTop]
#    type = DirichletBC
#    variable = disp_y
#    boundary = top
#    value = 0.0
#  []
  [fixBottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
    preset=true
  []
  [Pressure]
    [inside]
      boundary = left
      function = inner_pressure
    []
  []
[]

[Materials]

  [./stress]
    type = NEMLStress
    database = 'gr91_v2_stochastic.xml'
    model = 'vonmises'
    neml_variable_iname = 'var0 var1 var2 var3 var4 var5 var6 var7'
    neml_variable_value0=8.274
    neml_variable_value1=747.4
    neml_variable_value2=3.55
    neml_variable_value3=112.0
    neml_variable_value4=44.33
    neml_variable_value5=650.6
    neml_variable_value6=10.71
    neml_variable_value7=2.04
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  line_search = 'none'

  l_max_its = 60
  l_tol = 8e-3

  nl_max_its = 40
  nl_rel_tol = 5e-4
  nl_abs_tol = 1e-7

  end_time = 6.3e8
  dtmin = 1
  dtmax = 1e8
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 100
    iteration_window = 4
    optimal_iterations = 12
  [../]
  [./Predictor]
    type = SimplePredictor
    scale = 1.0
  [../]
[]

[Postprocessors]
  [./num_lin_it]
    type = NumLinearIterations
  [../]
  [./num_nonlin_it]
    type = NumNonlinearIterations
  [../]
  [./tot_lin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_lin_it
  [../]
  [./tot_nonlin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_nonlin_it
  [../]


  [./end_time]
    type = 'FunctionValuePostprocessor'
    function = timefunc
  [../]


  [./max_strain_xx]
    type = ElementExtremeValue
    variable = strain_xx
    value_type = max
  [../]
  [./max_strain_yy]
    type = ElementExtremeValue
    variable = strain_yy
    value_type = max
  [../]
  [./max_strain_zz]
    type = ElementExtremeValue
    variable = strain_zz
    value_type = max
  [../]

  [./max_principal_strain]
    type = ElementExtremeValue
    variable = max_principal_strain
    value_type = max
  [../]
  [./max_principal_stress]
    type = ElementExtremeValue
    variable = max_principal_stress
    value_type = max
  [../]
[]

[UserObjects]
  [max_strain]
    type = Terminator
    expression = 'max_principal_strain > 0.01'
  []
[]

[Outputs]
  exodus = off
  csv = false
[]
