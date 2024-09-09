[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 1
  []
[]

[NEML2]
  input = 'models/polynomial_rom3.i'
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
  [ep]
    initial_condition = 0
  []
  [s1]
    initial_condition = 0
  []
  [s2]
    initial_condition = 0
  []
  [s3]
    initial_condition = 0
  []

[]

# //"in_stress": [0,...,300]
# //"in_temperature": [600,...,1100]
# // "in_plastic_strain": [0,...,0.01]
# // "in_cell": [1.485051e6, 8.5e12]
# // "in_wall": [5.033e12, 1.2e13]
# // "in_env": [1e-9, 1e-6]

[AuxKernels]
  [stress]
    type = ParsedAux
    variable = s
    use_xyzt = true
    expression = '77*t'
  []
  [temperature]
    type = ParsedAux
    variable = T
    use_xyzt = true
    expression = '600+100*t'
  []
  [plastic_strain]
    type = ParsedAux
    variable = ep
    use_xyzt = true
    expression = '0.001*t'
  []
  [cell]
    type = ParsedAux
    variable = s1
    use_xyzt = true
    expression = '1e12*t'
  []
  [wall]
    type = ParsedAux
    variable = s2
    use_xyzt = true
    expression = '5e12+1e12*t'
  []
  [env]
    type = ParsedAux
    variable = s3
    use_xyzt = true
    expression = '1e-8*t'
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
  [ep]
    type = MOOSEVariableToNEML2
    moose_variable = ep
    neml2_variable = state/ep
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
  [s3]
    type = MOOSEVariableToNEML2
    moose_variable = s3
    neml2_variable = state/s3
  []

  [rom_model]
    type = ExecuteNEML2Model
    model = rom
    gather_uos = 'trial_stress temperature ep s1 s2 s3'
    enable_AD = true
  []
[]

[Materials]
  [neml2_ep_rate]
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

[Executioner]
  type = Transient
  num_steps = 2
  dt = 1
  dtmin = 1
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
