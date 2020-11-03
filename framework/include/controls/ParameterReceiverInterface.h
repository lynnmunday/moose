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
#include "MultiAppTransfer.h"

// Forward Declarations
class ParameterReceiver;
class FEProblemBase;

class ParameterReceiverInterface
{
public:
  /**
   * Constructor
   */
  ParameterReceiverInterface(const MooseObject * moose_object);

  /**
   * This class gives access to the transfer method for controllable values
   * @param receiver - ParameterReceiver for the transfer
   * @param names - controllable object names
   * @param values - value to give to controllable objects
   */
  void transferParameters(ParameterReceiver & receiver,
                          const std::vector<std::string> & names,
                          const std::vector<Real> & values) const;

  /**
   * Return the ParameterReceiver object and perform error checking.
   * @param app_index The global sup-app index
   */
  ParameterReceiver & getReceiver2(const std::string & multi_app_name,
                                   const std::string & receiver_name,
                                   unsigned int app_index);

private:
  /// Reference the FEProblemBase class
  FEProblemBase & _feproblem;
};
