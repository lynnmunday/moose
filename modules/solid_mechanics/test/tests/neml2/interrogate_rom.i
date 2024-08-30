[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
  []
[]

[NEML2]
  input = 'models/polynomial_rom.i'
  model = 'rom'
  verbose = true
  mode = PARSE_ONLY
  device = 'cpu'
  enable_AD = true
[]

[AuxVariables]
  [s]
    initial_condition = 0
  []
  [T]
    initial_condition = 0
  []
  [s1]
    initial_condition = 0
  []
  [s2]
    initial_condition = 0
  []
[]
[AuxKernels]
  [s]
    type = ParsedAux
    variable = s
    use_xyzt = true
    expression = '50*t'
  []
  [T]
    type = ParsedAux
    variable = T
    use_xyzt = true
    expression = '300*t'
  []
[]

[Variables]
  [u]
  []
[]

[Problem]
  kernel_coverage_check = false
[]

[UserObjects]
  [trial_stress]
    type = MOOSEVariableToNEML2
    moose_variable = s
    neml2_variable = state/s
  []
  [temperature]
    type = MOOSEVariableToNEML2
    moose_variable = T
    neml2_variable = forces/T
  []
  [s1]
    type = MOOSEVariableToNEML2
    moose_variable = s1
    neml2_variable = state/s1
  []
  [s2]
    type = MOOSEVariableToNEML2
    moose_variable = s2
    neml2_variable = state/s2
  []

  [rom_model]
    type = ExecuteNEML2Model
    model = rom
    gather_uos = 'trial_stress temperature s1 s2'
    enable_AD = true
  []
[]

[Materials]
  [neml2_3p_rate]
    type = NEML2ToRealMOOSEMaterialProperty
    execute_neml2_model_uo = rom_model
    neml2_variable = state/ep_rate
    moose_material_property = ep_rate
    outputs = exodus
  []
  [neml2_s1_rate]
    type = NEML2ToRealMOOSEMaterialProperty
    execute_neml2_model_uo = rom_model
    neml2_variable = state/s1_rate
    moose_material_property = s1_rate
    outputs = exodus
  []
  [neml2_s2_rate]
    type = NEML2ToRealMOOSEMaterialProperty
    execute_neml2_model_uo = rom_model
    neml2_variable = state/s2_rate
    moose_material_property = s2_rate
    outputs = exodus
  []
[]

[AuxVariables]
  [T]
  []
[]

[Executioner]
  type = Transient
  num_steps = 2
[]

[Postprocessors]
  [ep_rate]
    type = ElementAverageMaterialProperty
    mat_prop = ep_rate
  []
  [s1_rate]
    type = ElementAverageMaterialProperty
    mat_prop = s1_rate
  []
  [s2_rate]
    type = ElementAverageMaterialProperty
    mat_prop = s2_rate
  []
[]

[Outputs]
  exodus = true
[]
