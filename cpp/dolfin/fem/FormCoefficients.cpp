// Copyright (C) 2018 Garth N. Wells
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "FormCoefficients.h"
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/log/log.h>
#include <memory>
#include <string>
#include <ufc.h>

using namespace dolfin;
using namespace dolfin::fem;

//-----------------------------------------------------------------------------
FormCoefficients::FormCoefficients(const ufc_form& ufc_form)
    : _coefficients(ufc_form.num_coefficients)
{
  // Create finite elements for coefficients
  for (int i = 0; i < ufc_form.num_coefficients; i++)
    _original_pos.push_back(ufc_form.original_coefficient_position(i));
}
//-----------------------------------------------------------------------------
std::size_t FormCoefficients::size() const { return _coefficients.size(); }
//-----------------------------------------------------------------------------
std::vector<int> FormCoefficients::offsets() const
{
  std::vector<int> n = {0};
  for (auto& c : _coefficients)
  {
    if (!c)
      throw std::runtime_error("Not all form coefficients have been set.");
    n.push_back(n.back() + c->function_space()->element()->space_dimension());
  }
  return n;
}
//-----------------------------------------------------------------------------
void FormCoefficients::set(
    std::size_t i, std::shared_ptr<const function::Function> coefficient)
{
  assert(i < _coefficients.size());
  _coefficients[i] = coefficient;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const function::Function>
FormCoefficients::get(std::size_t i) const
{
  assert(i < _coefficients.size());
  return _coefficients[i];
}
//-----------------------------------------------------------------------------
std::size_t FormCoefficients::original_position(std::size_t i) const
{
  assert(i < _original_pos.size());
  return _original_pos[i];
}
//-----------------------------------------------------------------------------
