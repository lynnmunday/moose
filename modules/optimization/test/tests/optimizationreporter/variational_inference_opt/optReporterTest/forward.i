[Mesh]
  type = GeneratedMesh
  dim = 1
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[AuxVariables]
  [u]
  []
[]

[Reporters]
  [forward]
    type = ConstantReporter
    real_vector_names = 'p1 p2'
    real_vector_values = '0.1 0.2; 0.3 0.4 0.5'
    real_names = 'obj'
    real_values = '0.5'
  []
[]

[Outputs]
[]
