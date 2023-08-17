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
#include "AuxiliarySystem.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

registerMooseObject("OptimizationApp", NodalResidualDerivative);

InputParameters
NodalResidualDerivative::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription("Multiply off-diag Jacobian by variable");

  params.addRequiredParam<VariableName>("u_var",
                                        "Variable being differentiated (dependent variable)");
  params.addRequiredParam<VariableName>(
      "p_var", "Differentiated with respect to variable (independent variable)");
  params.addRequiredParam<AuxVariableName>("dudp_var", "Target aux");
  params.addRequiredParam<std::string>("matrix_tag", "Matrix Tag");
  params.addParam<NonlinearSystemName>(
      "nl_sys_name", "nl0", "Name of nonlinear system containing coupled variable");

  return params;
}

NodalResidualDerivative::NodalResidualDerivative(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _nl_sys(_fe_problem.getNonlinearSystemBase()),
    _aux_sys(_fe_problem.getAuxiliarySystem()),
    _nl_sys_num(_fe_problem.nlSysNum(getParam<NonlinearSystemName>("nl_sys_name"))),
    _aux_sys_num(_aux_sys.number()),
    _u_var(_subproblem.getStandardVariable(_tid, getParam<VariableName>("u_var"))),
    _p_var(_subproblem.getStandardVariable(_tid, getParam<VariableName>("p_var"))),
    _dudp_var(_subproblem.getStandardVariable(_tid, getParam<AuxVariableName>("dudp_var"))),
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
  const auto & mesh = _fe_problem.mesh();
  const auto & aslr = *mesh.getActiveSemiLocalNodeRange();
  std::vector<std::tuple<dof_id_type, dof_id_type, dof_id_type>> dofs(aslr.size());

  // collect all dofs
  const auto n_dofs = aslr.size();
  std::vector<dof_id_type> u_dofs(n_dofs), p_dofs(n_dofs), dudp_dofs(n_dofs);
  std::size_t i = 0;
  for (const auto node : aslr)
  {
    u_dofs[i] = node->dof_number(_nl_sys_num, _u_var.number(), 0);
    p_dofs[i] = node->dof_number(_nl_sys_num, _p_var.number(), 0);
    dudp_dofs[i] = node->dof_number(_aux_sys_num, _dudp_var.number(), 0);
    ++i;
  }

  // get u solution
  std::vector<Real> u(n_dofs);
  _nl_sys.currentSolution()->get(u_dofs, u);

  // get entire  jacobian
  SparseMatrix<Number> & jacobian = _nl_sys.getMatrix(_tag_id);

  // get R_up subblock
  PetscMatrix<Number> R_up(_communicator);
  jacobian.create_submatrix(R_up, u_dofs, p_dofs);

  std::cout << R_up << '\n';

  // PetscVector<Number> dest(_communicator);
  // R_up.vector_mult(dest, *_nl_sys.currentSolution());

  // std::cout << dest << std::endl;
  // // maybe do something like in PhysicsBasedPreconditioner::setup()
}
