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
  void fillInitialCovarianceParamsVector(std::vector<Real> & initial_conditions);
  /**
   * Function to fill the covariance vector of parameters bounds.
   */
  void fillBoundsCovarianceParamsVector(std::vector<Real> & data_vec, Real value);

private:
  /// seed for random number generator
  unsigned int _seed;
};
