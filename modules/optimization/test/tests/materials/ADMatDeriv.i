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
  [diff]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxVariables]
  [diff_star]
    family = MONOMIAL
    order = CONSTANT
  []
  [d_mat_d_diff]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [d_mat_d_diff]
      type = ADMaterialRealAux
      variable = d_mat_d_diff
      property = ad_opt_prop_deriv
      execute_on = 'TIMESTEP_END'
    []
[]

[ICs]
  [set_ic]
    type = FunctionIC
    variable = diff_star
    function = '3+x+y'
  []
[]

[Materials]
  [difference]
    type = ADParsedMaterial
    expression = 'diff_star'
    coupled_variables = 'diff_star'
    property_name = difference
  []
  [ADmaterial_property]
    type = ADParsedMaterial
    expression = "diff * diff"
    coupled_variables = "diff"
    property_name = ad_opt_prop
  []
  [ADmaterialDeriv_property]
    type = ADCoupledVariableDerivativeMaterial
    derivative_prop_name = ad_opt_prop_deriv
    material_property = ad_opt_prop
    derivative_variable_name = diff
  []
[]

[Kernels]
  [diffusion]
    type = ADMaterialPropertyValue
    variable = diff
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
