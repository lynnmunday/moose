//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "ParameterReceiverInterface.h"
#include "ParameterReceiver.h"
#include "FEProblemBase.h"
#include "MultiApp.h"

// could initilaize list of multiapp references.
ParameterReceiverInterface::ParameterReceiverInterface(const MooseObject * moose_object)
  : _feproblem(
        *(moose_object->parameters()).getCheckedPointerParam<FEProblemBase *>("_fe_problem_base"))
{
}

void
ParameterReceiverInterface::transferParameters(ParameterReceiver & receiver,
                                               const std::vector<std::string> & names,
                                               const std::vector<Real> & values) const
{
  receiver.transfer(names, values);
}

ParameterReceiver &
ParameterReceiverInterface::getReceiver2(const std::string & multi_app_name,
                                         const std::string & receiver_name,
                                         unsigned int app_index)
{
  // Test that the sub-application has the given Control object
  std::shared_ptr<MultiApp> multi_app = _feproblem.getMultiApp(multi_app_name);
  FEProblemBase & to_problem = multi_app->appProblemBase(app_index);
  ExecuteMooseObjectWarehouse<Control> & control_wh = to_problem.getControlWarehouse();
  if (!control_wh.hasActiveObject(receiver_name))
    mooseError("The sub-application (",
               multi_app->name(),
               ") does not contain a Control object with the name '",
               receiver_name,
               "'.");

  ParameterReceiver * ptr =
      dynamic_cast<ParameterReceiver *>(control_wh.getActiveObject(receiver_name).get());

  if (!ptr)
    mooseError(
        "The sub-application (",
        multi_app->name(),
        ") Control object for the 'to_control' parameter must be of type 'ParameterReceiver'.");

  return *ptr;
}
