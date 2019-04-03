// Copyright (C) 2019 Igor A. Baratta
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "PointSource.h"
#include <dolfin/geometry/Point.h>
#include <dolfin/function/FunctionSpace.h>


using namespace dolfin;



//-----------------------------------------------------------------------------
fem::PointSource::PointSource(std::shared_ptr<const function::FunctionSpace> V, 
                              const geometry::Point& point,
                                PetscScalar magnitude)
{


}
//-----------------------------------------------------------------------------