//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADCoupledVariableDerivativeMaterial.h"

#include "metaphysicl/raw_type.h"
#include "NonlinearSystem.h"

registerMooseObject("OptimizationApp", ADCoupledVariableDerivativeMaterial);

InputParameters
ADCoupledVariableDerivativeMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>(
      "derivative_prop_name",
      "The name of the material property where we store the derivative values.");
  params.addRequiredParam<MaterialPropertyName>(
      "material_property", "Name of material property being differentiated by coupled_variable.");
  params.addRequiredParam<std::string>(
      "derivative_variable_name",
      "The variable name that material_property will be differentiated by.");
  params.addParam<NonlinearSystemName>(
      "nl_sys_name", "nl0", "Name of nonlinear system containing coupled variable");
  params.addClassDescription("Stores derivative of a variable into material properties");
  return params;
}

ADCoupledVariableDerivativeMaterial::ADCoupledVariableDerivativeMaterial(
    const InputParameters & parameters)
  : Material(parameters),
    _nl_sys(_fe_problem.getNonlinearSystemBase()),
    _nl_sys_num(_fe_problem.nlSysNum(getParam<NonlinearSystemName>("nl_sys_name"))),
    _dmat_dparam(declareADProperty<Real>("derivative_prop_name")),
    _mat(getADMaterialProperty<Real>("material_property")),
    _var_num(_nl_sys.getVariable(_nl_sys_num, getParam<std::string>("derivative_variable_name"))
                 .number())
{
}

void
ADCoupledVariableDerivativeMaterial::computeQpProperties()
{
  dof_id_type dof_id = _current_elem->dof_number(_nl_sys_num, _var_num, /*tid*/ 0);
  _dmat_dparam[_qp] = _mat[_qp].derivatives()[dof_id];
}
