// Copyright (C) 2019 Igor A. Baratta
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#pragma once

#include <dolfin/geometry/Point.h>
#include <petscsys.h>
#include <memory>

namespace dolfin
{

namespace function
{
class FunctionSpace;
}

namespace la
{
class PetscVector;
}


namespace fem
{

class PointSource
{
    public:
    /// Create point source at given point of given magnitude
    PointSource(std::shared_ptr<const function::FunctionSpace> V,
                const geometry::Point& point,
                PetscScalar magnitude);
    


};

}

}