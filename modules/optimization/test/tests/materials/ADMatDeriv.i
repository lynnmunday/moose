[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
    elem_type = QUAD4
  []
[]

[Variables]
  [coeff]
    family = MONOMIAL
    order = CONSTANT
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
    property = ad_opt_prop_deriv
    execute_on = 'TIMESTEP_END'
  []
[]

[ICs]
  # this is the parameter coeff that TAO is going to set
  [coeff_from_transfer]
    type = FunctionIC
    variable = controllable_coeff
    function = '3+x+y'
  []
[]

[Materials]
  [difference]
    type = ADParsedMaterial
    expression = 'controllable_coeff'
    coupled_variables = 'controllable_coeff'
    property_name = difference
  []
  [ADmaterial_property]
    type = ADParsedMaterial
    expression = "coeff * coeff"
    coupled_variables = "coeff"
    property_name = ad_opt_prop
  []
  [ADmaterialDeriv_property]
    type = ADCoupledVariableDerivativeMaterial
    derivative_prop_name = ad_opt_prop_deriv
    material_property = ad_opt_prop
    derivative_variable_name = coeff
  []
[]

[Kernels]
  [diffusion]
    type = ADMaterialPropertyValue
    variable = coeff
    prop_name = difference
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = none
[]

[Outputs]
  exodus = true
[]
