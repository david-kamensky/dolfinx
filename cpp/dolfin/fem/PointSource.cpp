// Copyright (C) 2019 Igor A. Baratta
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "PointSource.h"
#include <dolfin/geometry/Point.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/PETScVector.h>


using namespace dolfin;



//-----------------------------------------------------------------------------
fem::PointSource::PointSource(std::shared_ptr<function::FunctionSpace> V,
                              const geometry::Point& point,
                              PetscScalar magnitude): _V0(V)
{
    // Check if FunctionSpace is scalar
    // FIXME: Allow vectoria Point Sources
    assert(V->element()->space_dimension() == 1);
    _point = point;
    _magniude = magnitude;
}
//-----------------------------------------------------------------------------

void fem::PointSource::apply(la::PETScVector& b)
{
    assert(b.size() > 0);

}