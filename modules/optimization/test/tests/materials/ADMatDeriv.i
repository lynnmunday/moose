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
  #dummy variable of parameter
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
    property = mat_prop_deriv
    execute_on = 'TIMESTEP_END'
  []
[]

[ICs]
  # this is the parameter coeff that TAO is going to set
  [coeff_from_transfer]
    type = FunctionIC
    variable = controllable_coeff
    function = 'x+y'
  []
[]

[Materials]
  [controllable_coeff_material]
    #This puts the controllable_coeff AuxVariable set by an IC and places it in a material
    type = ADParsedMaterial
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
    type = ADMaterialPropertyValue
    variable = coeff
    prop_name = controllable_coeff_material
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
