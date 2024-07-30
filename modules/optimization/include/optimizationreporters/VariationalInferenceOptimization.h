//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralOptimization.h"

/**
 * Optimization reporter that interfaces with TAO. Objective value, gradients,
 * and constraints need to calculated in the subapps.
 */
class VariationalInferenceOptimization : public GeneralOptimization
{
public:
  static InputParameters validParams();
  VariationalInferenceOptimization(const InputParameters & parameters);

  virtual Real computeObjective() override;

protected:
  /// sampled parameter values declared as reporter data used in the forward problem
  std::vector<std::vector<Real> *> _sampled_parameters;

  virtual void setICsandBounds() override;
  /**
   * Function to set parameters.
   * This is the first function called in objective/gradient/hessian routine
   */
  virtual void updateParameters(const libMesh::PetscVector<Number> & x) override;
  /**
   * Function to fill the covariance vector of parameters initial conditions.
   */
  std::vector<Real> fillInitialCovarianceParamsVector();
  /**
   * Function to fill the covariance vector of parameters bounds.
   */
  void fillBoundsCovarianceParamsVector(std::vector<Real> & data_vec, Real value);

private:
  /// Reporter that will hold the ELBO objective value
  Real & _elbo_objective_val;
  /// seed for random number generator
  const unsigned int _seed;
  const Real _experimental_noise;
  const unsigned int _num_experiments;
  const unsigned int _num_measurements_per_experiment;
  std::vector<std::vector<Real>> _random_draws;
  Real _const_term_loglikelihood;
};
