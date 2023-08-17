[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
    xmin = 0
    xmax = 2
    ymin = 0
    ymax = 2
    elem_type = QUAD4
  []
  [nodeset_all]
    type = BoundingBoxNodeSetGenerator
    input = gmg
    bottom_left = '-0.1 -0.1 0'
    top_right = '2.1 2.1 0'
    new_boundary = 'all_nodes'
  []
[]

[Problem]
  kernel_coverage_check = false
  extra_tag_matrices = 'mat_tag'
[]

[UserObjects]
  [nodalDer]
    type = NodalResidualDerivative
    u_var = forwardT
    p_var = coeff
    dudp_var = dudp
    execute_on = TIMESTEP_END
    matrix_tag = mat_tag
  []
[]

[Variables]
  [forwardT]
  []
[]
[AuxVariables]
  [dudp][]
[]
[Kernels]
  [heat_conduction]
    type = ADMatDiffusion
    variable = forwardT
    diffusivity = 'mat_prop'
    extra_matrix_tags = 'mat_tag'
  []
  [heat_source]
    type = ADBodyForce
    function = volumetric_heat_func
    variable = forwardT
    extra_matrix_tags = 'mat_tag'
  []
[]

[Reporters]
  [params]
    type = ConstantReporter
    real_vector_names = 'heat_source'
    real_vector_values = '333' # Dummy
  []
[]
[Functions]
  [volumetric_heat_func]
    type = ParsedOptimizationFunction
    expression = q
    param_symbol_names = 'q'
    param_vector_name = 'params/heat_source'
  []
[]
[BCs]
  [left]
    type = NeumannBC
    variable = forwardT
    boundary = left
    value = 0
    extra_matrix_tags = 'mat_tag'
  []
  [right]
    type = NeumannBC
    variable = forwardT
    boundary = right
    value = 0
    extra_matrix_tags = 'mat_tag'
  []
  [bottom]
    type = DirichletBC
    variable = forwardT
    boundary = bottom
    value = 2
    extra_matrix_tags = 'mat_tag'
  []
  [top]
    type = DirichletBC
    variable = forwardT
    boundary = top
    value = 1
    extra_matrix_tags = 'mat_tag'
  []
[]

[Preconditioning]
  active = 'FSP'
  [FSP]
    type = FSP
    full = false
    topsplit = 'forwardT_coeff'
    [forwardT_coeff]
      splitting = 'forwardT coeff'
      # Generally speaking, there are four types of splitting we could choose
      # <additive,multiplicative,symmetric_multiplicative,schur>
      splitting_type = additive
    []
    [forwardT]
      vars = 'forwardT'
      # PETSc options for this subsolver
      # A prefix will be applied, so just put the options for this subsolver only
      petsc_options_iname = '-ksp_type -pc_type -pc_factor_shift_type'
      petsc_options_value = 'preonly lu nonzero'
    []
    [coeff]
      vars = 'coeff'
      # PETSc options for this subsolver
      petsc_options_iname = '-ksp_type -pc_type -pc_factor_shift_type'
      petsc_options_value = 'preonly lu nonzero'
    []
  []
[SMP]
  type = SMP
  full = true
  petsc_options = '-ksp_view_pmat'
[]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  # automatic_scaling = true  #check that this doesn't change the answer.... i think it will because it changes the jacobian
  line_search = none
  nl_forced_its = 1
  petsc_options = '-ksp_view_pmat'
[]

[Outputs]
  exodus = true
[]

#-------------- material derivative ------------------------

[Variables]
  #dummy variable of parameter
  [coeff]
    family = LAGRANGE
    order = FIRST
    initial_condition = 1
  []
[]

[BCs]
  [set_coeff]
    type = DirichletBC
    variable = coeff
    value = '2'
    boundary = all_nodes
  []
[]

[Materials]
  [ADmaterial_property]
    # this is just a new type material that is parameter^2 so that we have something to differentiate
    # This is the material property being used
    type = ADParsedMaterial
    expression = "coeff * coeff"
    coupled_variables = "coeff"
    property_name = mat_prop
  []
[]

