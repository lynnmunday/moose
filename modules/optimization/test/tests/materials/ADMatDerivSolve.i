[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
    xmin = 0
    xmax = .6
    ymin = 0
    ymax = 2
    elem_type = QUAD4
  []
[]

[Variables]
  [forwardT]
  []
[]
[Kernels]
  [heat_conduction]
    type = ADMatDiffusion
    variable = forwardT
    diffusivity = 'mat_prop'
  []
  [heat_source]
    type = ADBodyForce
    function = volumetric_heat_func
    variable = forwardT
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
  []
  [right]
    type = NeumannBC
    variable = forwardT
    boundary = right
    value = 0
  []
  [bottom]
    type = DirichletBC
    variable = forwardT
    boundary = bottom
    value = 2
  []
  [top]
    type = DirichletBC
    variable = forwardT
    boundary = top
    value = 1
  []
[]



#-------------- material derivative ------------------------
[Variables]
  #dummy variable of parameter
  [coeff]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 1
  []
[]

[AuxVariables]
  [controllable_coeff]
    family = MONOMIAL
    order = CONSTANT
  []
  [d_mat_d_coeff]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [d_mat_d_coeff]
    type = ADMaterialRealAux
    variable = d_mat_d_coeff
    property = mat_prop_deriv
    execute_on = 'TIMESTEP_END'
  []
[]

[ICs]
  # this is the parameter coeff that TAO is going to set
  [coeff_from_transfer]
    type = FunctionIC
    variable = controllable_coeff
    function = '3'
  []
[]

[Materials]
  [controllable_coeff_material]
    #This puts the controllable_coeff AuxVariable set by an IC and places it in a material
    type = ParsedMaterial
    expression = 'controllable_coeff'
    coupled_variables = 'controllable_coeff'
    property_name = controllable_coeff_material
  []
  [ADmaterial_property]
    # this is just a new type material that is parameter^2 so that we have something to differentiate
    # This is the material property being used
    type = ADParsedMaterial
    expression = "coeff * coeff"
    coupled_variables = "coeff"
    property_name = mat_prop
  []
  [ADmaterialDeriv_property]
    type = ADCoupledVariableDerivativeMaterial
    derivative_prop_name = mat_prop_deriv
    material_property = mat_prop
    derivative_variable_name = coeff
  []
[]

[Kernels]
  [setCoefficient]
    type = MaterialPropertyValue
    variable = coeff
    prop_name = controllable_coeff_material
  []
[]

[Preconditioning]
    active = 'SMP'
    # [FSP]
    #   type = FSP
    #   full = false
    #   topsplit = 'forwardT_coeff'
    #   [forwardT_coeff]
    #     splitting = 'forwardT coeff'
    #     # Generally speaking, there are four types of splitting we could choose
    #     # <additive,multiplicative,symmetric_multiplicative,schur>
    #     splitting_type = additive
    #   []
    #   [forwardT]
    #     vars = 'forwardT'
    #     # PETSc options for this subsolver
    #     # A prefix will be applied, so just put the options for this subsolver only
    #     petsc_options_iname = '-ksp_type -pc_type -pc_factor_shift_type'
    #     petsc_options_value = 'preonly lu nonzero'
    #   []
    #   [coeff]
    #     vars = 'coeff'
    #     # PETSc options for this subsolver
    #     petsc_options_iname = '-ksp_type -pc_type -pc_factor_shift_type'
    #     petsc_options_value = 'preonly lu nonzero'
    #   []
    # []
  [SMP]
    type = SMP
    full = true
    petsc_options = '-ksp_view_pmat'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  automatic_scaling = true
  line_search = none
  nl_forced_its = 2
[]

[Outputs]
  exodus = true
[]
