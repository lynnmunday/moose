# To run this input file with a normal tensor mechanics executable, use the command
#  tensor_mechanics-opt -i fail_negative_trial_stress_2drz.i --allow-test-objects

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

[AuxVariables]
  [temperature]
    initial_condition = 900.0
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
        use_automatic_differentiation = true
        strain = FINITE
        add_variables = true
        generate_output = 'max_principal_stress max_principal_strain vonmises_stress
                          strain_xx strain_yy strain_zz
                          stress_xx stress_yy stress_zz'
        extra_vector_tags = 'ref'
      []
    []
  []
[]

[Functions]
  [timefunc]
    type=ParsedFunction
    value = t
  []
[]

[BCs]
#  [fixTop]
#    type = DirichletBC
#    variable = disp_y
#    boundary = top
#    value = 0.0
#  []
  [fixBottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
    preset=true
  []
  [InsidePressure]
    type = ADPressure
    boundary = left
    variable = disp_x
    component = 0
    constant = 2.5e6
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 3.30e11
    poissons_ratio = 0.3
  []
  [./stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = rom_stress_prediction
  [../]
  [./rom_stress_prediction]
#    internal_solve_output_on = always
    type = SS316HLAROMANCEStressUpdateTest
    temperature = temperature
    extrapolate_to_zero_stress = true
    initial_mobile_dislocation_density = 9.0e12
    initial_immobile_dislocation_density = 8.6e11
    outputs = all
    input_window_failure_action = ignore
#    verbose=true
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
  [./effective_strain_avg]
    type = ElementAverageValue
    variable = effective_creep_strain
  [../]
  [./temperature]
    type = ElementAverageValue
    variable = temperature
  [../]
  [./mobile_dislocations]
    type = ElementAverageValue
    variable = mobile_dislocations
  [../]
  [./immobile_disloactions]
    type = ElementAverageValue
    variable = immobile_dislocations
  [../]

  [./num_lin_it]
    type = NumLinearIterations
  [../]
  [./num_nonlin_it]
    type = NumNonlinearIterations
  [../]

  [./end_time]
    type = 'FunctionValuePostprocessor'
    function = timefunc
  [../]

  [./max_vonmises_stress]
    type = ElementExtremeValue
    variable = vonmises_stress
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
