//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NodalResidualDerivative.h"
#include "NonlinearSystemBase.h"

registerMooseObject("OptimizationApp", NodalResidualDerivative);

InputParameters
NodalResidualDerivative::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription("Multiply off-diag Jacobian by variable");

  params.addRequiredParam<VariableName>("var",
                                        "Variable being differentiated (dependent variable)");
  params.addRequiredParam<VariableName>(
      "deriv_wrt_var", "Differentiated with respect to variable (independent variable)");
  params.addRequiredParam<std::string>("matrix_tag", "Matrix Tag");
  params.addParam<NonlinearSystemName>(
      "nl_sys_name", "nl0", "Name of nonlinear system containing coupled variable");

  return params;
}

NodalResidualDerivative::NodalResidualDerivative(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _nl_sys(_fe_problem.getNonlinearSystemBase()),
    _nl_sys_num(_fe_problem.nlSysNum(getParam<NonlinearSystemName>("nl_sys_name"))),
    _var(_subproblem.getStandardVariable(_tid, getParam<VariableName>("var"))),
    _deriv_wrt_var(_subproblem.getStandardVariable(_tid, getParam<VariableName>("deriv_wrt_var"))),
    _tag_id(_subproblem.getMatrixTagID(getParam<std::string>("matrix_tag")))
{
}

void
NodalResidualDerivative::initialize()
{
}

void
NodalResidualDerivative::execute()
{
  dof_id_type dof_id_var = _current_node->dof_number(_nl_sys_num, _var.number(), _tid);
  dof_id_type dof_id_dvVar = _current_node->dof_number(_nl_sys_num, _deriv_wrt_var.number(), _tid);

  SparseMatrix<Number> & jacobian = _nl_sys.getMatrix(_tag_id); //_nl_sys.systemMatrixTag());

  std::cout << "dof_id_var,dof_id_dvVar= " << dof_id_var << " " << dof_id_dvVar << " "
            << jacobian(dof_id_var, dof_id_dvVar) << std::endl;
  // maybe do something like in PhysicsBasedPreconditioner::setup()
}
