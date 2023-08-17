//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "GeneralUserObject.h"

/**
 * User object that multiplies variable A by the off-diagaonal Jacobian term of variable A with
 * respect to variable B.
 */
class NodalResidualDerivative : public GeneralUserObject
{
public:
  static InputParameters validParams();
  NodalResidualDerivative(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override {}

private:
  NonlinearSystemBase & _nl_sys;
  AuxiliarySystem & _aux_sys;
  /// The number of the nonlinear system containing the coupled variable
  const unsigned int _nl_sys_num;
  /// The number of the nonlinear system containing the coupled variable
  const unsigned int _aux_sys_num;

  // Variable being differentiated (dependent variable)
  const MooseVariable & _u_var;
  // differentiated with respect to variable (independent variable)
  const MooseVariable & _p_var;

  // target variable
  const MooseVariable & _dudp_var;

  const unsigned int _tag_id;
};
