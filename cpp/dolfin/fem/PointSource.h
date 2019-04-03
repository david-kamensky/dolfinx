// Copyright (C) 2019 Igor A. Baratta
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#pragma once

#include <dolfin/geometry/Point.h>
#include <memory>
#include <petscsys.h>
#include <utility>

namespace dolfin
{

namespace function
{
class FunctionSpace;
}

namespace la
{
class PETScVector;
}

namespace fem
{

class PointSource
{
public:
  /// Create point source at given point of given magnitude
  PointSource(std::shared_ptr<function::FunctionSpace> V,
              const geometry::Point& point, PetscScalar magnitude);

  /// Apply point source to right-hand side vector
  void apply(la::PETScVector& b);

  /// Destructor
  ~PointSource() = default;

private:
  // Source term - pair of points and magnitude
  geometry::Point _point;
  PetscScalar _magnitude;
  std::shared_ptr<function::FunctionSpace> _V0;
};

} // namespace fem
} // namespace dolfin
