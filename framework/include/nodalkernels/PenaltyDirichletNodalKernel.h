//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NodalKernel.h"

// Forward Declarations
class PenaltyDirichletNodalKernel;

template <>
InputParameters validParams<PenaltyDirichletNodalKernel>();

class PenaltyDirichletNodalKernel : public NodalKernel
{
public:
  PenaltyDirichletNodalKernel(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:
  const Real & _p;
  const Real & _v;
};