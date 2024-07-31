[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Postprocessors]
  [max]
    type = ElementExtremeValue
    variable = u
  []
  [min]
    type = ElementExtremeValue
    variable = u
    value_type = min
  []
[]
[Reporters]
  [reporter]
    type = ConstantReporter
    real_names = 'a b'
    real_values = '1 10'
    outputs = none
  []
[]

[VectorPostprocessors]
  [min_max]
    type = VectorOfPostprocessors
    postprocessors = 'min max'
    reporters = 'reporter/a reporter/b'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = false
  csv = true
[]
