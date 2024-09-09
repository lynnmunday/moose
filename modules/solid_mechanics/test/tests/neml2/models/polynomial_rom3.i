[Models]
  [rom]
    type = TabulatedPolynomialModel3
    model_file_name = '/Users/mundlb/projects/bison/data/laromance/lafleur_ht9_6Dgrid.json'
    von_mises_stress = 'state/s'
    temperature = 'forces/T'
    equivalent_plastic_strain = 'state/ep'
    internal_state_1 = 'state/s1'
    internal_state_2 = 'state/s2'
    internal_state_3 = 'state/s3'
    equivalent_plastic_strain_rate = 'state/ep_rate'
    internal_state_1_rate = 'state/s1_rate'
    internal_state_2_rate = 'state/s2_rate'

    # use_AD_first_derivative = true
    # use_AD_second_derivative = true
  []
[]
