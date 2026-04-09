[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 8
  []
[]

[AuxVariables]
  [aux]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [threaded_param_error]
    type = ThreadedParamErrorAux
    variable = aux
    trigger = true
    execute_on = INITIAL
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]
