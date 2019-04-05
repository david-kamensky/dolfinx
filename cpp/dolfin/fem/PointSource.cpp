// Copyright (C) 2019 Igor A. Baratta
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "PointSource.h"
#include <dolfin/common/types.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/geometry/BoundingBoxTree.h>
#include <dolfin/geometry/Point.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/mesh/Mesh.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
fem::PointSource::PointSource(std::shared_ptr<function::FunctionSpace> V,
                              const geometry::Point& point,
                              PetscScalar magnitude)
    : _V0(V)
{
  assert(V);
  assert(magnitude);

  _point = point;
  _magnitude = magnitude;
}
//-----------------------------------------------------------------------------

void fem::PointSource::apply(la::PETScVector& b)
{
  assert(b.size() > 0);

  const fem::FiniteElement& element = *_V0->element();

  // Get index of first cell containing the point
  const mesh::Mesh& mesh = *_V0->mesh();
  geometry::BoundingBoxTree tree(mesh, mesh.topology().dim());
  unsigned int cell_index = tree.compute_first_entity_collision(_point, mesh);

  // Create cell that contains the point
  mesh::Cell cell(mesh, cell_index);

  // Number of dofs per cell for the finite element
  std::size_t space_dimension = element.space_dimension();
  // Rank of the value space
  const std::size_t value_rank = element.value_rank();
  // Reference cell value size (1 for scalar)
  std::size_t reference_value_size = element.reference_value_size();

  // get cell dof coordinates (re-allocated inside function for thread safety)
  EigenRowArrayXXd coordinate_dofs(space_dimension, mesh.geometry().dim());
  cell.get_coordinate_dofs(coordinate_dofs);

  // Basis size (all dimensions)
  std::size_t size_basis = 1;
  for (std::size_t i = 0; i < value_rank; ++i)
  {
    // dimension of the value space for axis i
    size_basis *= element.value_dimension(i);
  }

  // Get cell dofmap
  const fem::GenericDofMap& dofmap = *_V0->dofmap();
  const Eigen::Map<const Eigen::Array<PetscInt, Eigen::Dynamic, 1>> cell_dofs
      = dofmap.cell_dofs(cell.index());

  // Data structure for basis functions at given point in reference cell
  // reference_values[num_points][space_demension][reference_value_size]
  Eigen::Tensor<double, 3, Eigen::RowMajor> reference_values(
      1, space_dimension, reference_value_size);

  std::vector<PetscScalar> values(space_dimension);
  std::vector<double> basis(size_basis);

  // FIXME: From geometry:Point to Eigen::VectorXd
  Eigen::VectorXd X(mesh.geometry().dim());
  for (std::size_t i = 0; i < mesh.geometry().dim(); i++)
    X(i) = _point.coordinates()[i];

  element.evaluate_reference_basis(reference_values, X);

  std::cout << "Imprime valores" << reference_values << '\n';

  for (std::size_t i = 0; i < space_dimension; ++i)
  {
    double basis_sum = 0.0;
    for (std::size_t j = 0; j < reference_value_size; ++j)
      basis_sum += reference_values(0, i, j);
    values[i] = _magnitude * basis_sum;
  }

  la::VecWrapper v(b.vec());
  for (std::size_t i = 0; i < space_dimension; ++i)
    v.x[cell_dofs[i]] = values[i];
  v.restore();
}
