//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ThreadedParamErrorAux.h"

registerMooseObject("MooseTestApp", ThreadedParamErrorAux);

InputParameters
ThreadedParamErrorAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addParam<bool>("trigger", true, "Parameter used to anchor the test paramError location.");
  return params;
}

ThreadedParamErrorAux::ThreadedParamErrorAux(const InputParameters & parameters)
  : AuxKernel(parameters)
{
  if (isNodal())
    paramError("variable", "This AuxKernel only supports Elemental fields");
}

Real
ThreadedParamErrorAux::computeValue()
{
  if (_current_elem->id() == 0)
    paramError("trigger", "Intentional threaded paramError from computeValue");

  return 0.0;
}
