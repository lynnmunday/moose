//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseTypes.h"
#include "NodalResidualDerivative.h"
#include "NonlinearSystemBase.h"
#include "AuxiliarySystem.h"
#include "Conversion.h"
#include "libmesh/id_types.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/sparse_matrix.h"
#include <algorithm>
#include <iostream>
#include <memory>

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
  // Now we need to get the row and column indicies entries we care about.
  std::vector<dof_id_type> row_idx;
  std::vector<dof_id_type> col_idx;

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

  // dudp vector (use aux var for this!)
  std::vector<Real> dudp_vec(n_dofs);

  // inverse p_dofs
  std::map<dof_id_type, std::size_t> inverse_p_dofs;
  for (const auto i : index_range(p_dofs))
    inverse_p_dofs[p_dofs[i]] = i;

  // get u solution
  std::vector<Real> us(n_dofs);
  _nl_sys.currentSolution()->get(u_dofs, us);

  // get entire  jacobian
  SparseMatrix<Number> & jacobian = _nl_sys.getMatrix(_tag_id);

  std::vector<Real> row;
  std::vector<dof_id_type> row_dofs;
  // u_idx is the row
  for (const auto u_idx : index_range(u_dofs))
  {
    const auto u_dof = u_dofs[u_idx];
    const auto u = us[u_idx];

    jacobian.get_row(u_dof, row_dofs, row);
    // j is the column?
    for (const auto j : index_range(row_dofs))
    {
      const auto it = inverse_p_dofs.find(row_dofs[j]);
      if (it == inverse_p_dofs.end())
        continue;

      const auto p_idx = it->second;
      const auto dudp = row[j];
      dudp_vec[p_idx] += dudp * u;
      // save the index information
      row_idx.push_back(u_idx);
      col_idx.push_back(j);
    }
  }

  // Here is the sparse matrix that we want
  libMesh::PetscMatrix<Real> dRdp(_communicator);
  // This will need to be different in parallel
  dRdp.init(u_dofs.size(), p_dofs.size(), u_dofs.size(), p_dofs.size());

  // Extract the sub components into the matrix
  dof_id_type idx = 0;
  for (MooseIndex(u_dofs.size()) r_idx = 0; r_idx < u_dofs.size(); r_idx++)
  {
    for (MooseIndex(p_dofs.size()) c_idx = 0; c_idx < p_dofs.size(); c_idx++)
    {
      dRdp.set(r_idx, c_idx, jacobian(row_idx[idx], col_idx[idx]));
      idx++;
    }
  }
  dRdp.close();

  for (size_t row = 0; row < dRdp.m(); row++)
  {
    std::vector<dof_id_type> row_ind;
    std::vector<Real> row_vals;
    dRdp.get_row(row, row_ind, row_vals);
    std::cout << Moose::stringify(row_vals) << '\n';
  }
  std::cout << Moose::stringify(dudp_vec) << '\n';
}
