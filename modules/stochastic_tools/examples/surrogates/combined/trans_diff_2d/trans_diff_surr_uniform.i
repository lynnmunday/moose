[StochasticTools]
[]

[Distributions]
  [C_dist]
    type = Uniform
    lower_bound = 0.01
    upper_bound = 0.02
  []
  [f_dist]
    type = Uniform
    lower_bound = 15
    upper_bound = 25
  []
  [init_dist]
    type = Uniform
    lower_bound = 270
    upper_bound = 330
  []
[]

[Samplers]
  [sample]
    type = LatinHypercube
    num_rows = 100000
    num_bins = 50
    distributions = 'C_dist f_dist init_dist'
    execute_on = PRE_MULTIAPP_SETUP
  []
[]

[Surrogates]
  [pc_max]
    type = PolynomialChaos
    filename = 'trainer_out_pc_max.rd'
  []
  [pc_min]
    type = PolynomialChaos
    filename = 'trainer_out_pc_min.rd'
  []
  [pr_min]
    type = PolynomialRegressionSurrogate
    filename = 'trainer_out_pr_min.rd'
  []
  [pr_max]
    type = PolynomialRegressionSurrogate
    filename = 'trainer_out_pr_max.rd'
  []
  [np_min]
    type = NearestPointSurrogate
    filename = 'trainer_out_np_min.rd'
  []
  [np_max]
    type = NearestPointSurrogate
    filename = 'trainer_out_np_max.rd'
  []
[]

# Computing statistics
[VectorPostprocessors]
  [pc_max_res]
    type = EvaluateSurrogate
    model = pc_max
    sampler = sample
  []
  [pc_max_stats]
    type = PolynomialChaosStatistics
    pc_name = 'pc_max'
    compute = 'mean stddev'
  []
  [pc_min_res]
    type = EvaluateSurrogate
    model = pc_min
    sampler = sample
  []
  [pc_min_stats]
    type = PolynomialChaosStatistics
    pc_name = 'pc_min'
    compute = 'mean stddev'
  []
  [pr_max_res]
    type = EvaluateSurrogate
    model = pr_max
    sampler = sample
  []
  [pr_max_stats]
    type = Statistics
    vectorpostprocessors = pr_max_res
    compute = 'mean stddev'
  []
  [pr_min_res]
    type = EvaluateSurrogate
    model = pr_min
    sampler = sample
  []
  [pr_min_stats]
    type = Statistics
    vectorpostprocessors = pr_min_res
    compute = 'mean stddev'
  []
  [np_max_res]
    type = EvaluateSurrogate
    model = np_max
    sampler = sample
  []
  [np_max_stats]
    type = Statistics
    vectorpostprocessors = np_max_res
    compute = 'mean stddev'
  []
  [np_min_res]
    type = EvaluateSurrogate
    model = np_min
    sampler = sample
  []
  [np_min_stats]
    type = Statistics
    vectorpostprocessors = np_min_res
    compute = 'mean stddev'
  []
[]

[Outputs]
  csv = true
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final'  # Default is "final"
    level = 4                    # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  []
[]