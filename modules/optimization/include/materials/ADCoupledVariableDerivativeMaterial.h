//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
class NonlinearSystemBase;
/**
 * Stores values of a variable into material properties
 */
class ADCoupledVariableDerivativeMaterial : public Material
{
public:
  static InputParameters validParams();

  ADCoupledVariableDerivativeMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  NonlinearSystemBase & _nl_sys;
  /// The number of the nonlinear system containing the coupled variable
  const unsigned _nl_sys_num;
  /// The name of the material property where the derivatives will be stored
  ADMaterialProperty<Real> & _dmat_dparam;
  /// material being differentiated by _param
  const ADMaterialProperty<Real> & _mat;
  /// coupled variable id for indexing into AD derivative
  const unsigned int _var_num;
};
