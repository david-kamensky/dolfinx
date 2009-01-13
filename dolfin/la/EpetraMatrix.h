// Copyright (C) 2008 Martin Sandve Alnes, Kent-Andre Mardal and Johannes Ring.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Garth N. Wells, 2008.
//
// First added:  2008-04-21
// Last changed: 2008-13-30

#ifndef __EPETRA_MATRIX_H
#define __EPETRA_MATRIX_H

#ifdef HAS_TRILINOS

#include <boost/shared_ptr.hpp>

#include <dolfin/log/LogStream.h>
#include <dolfin/common/Variable.h>
#include "GenericMatrix.h"

class Epetra_FECrsMatrix;
class Epetra_CrsGraph;

namespace dolfin
{

  class GenericSparsityPattern;

  /// This class provides a simple matrix class based on Epetra.
  /// It is a simple wrapper for an Epetra matrix object (Epetra_FECrsMatrix)
  /// implementing the GenericMatrix interface.
  ///
  /// The interface is intentionally simple. For advanced usage,
  /// access the Epetra_FECrsMatrix object using the function mat() and
  /// use the standard Epetra interface.

  class EpetraMatrix: public GenericMatrix, public Variable
  {
  public:

    /// Create empty matrix
    EpetraMatrix();

    /// Create M x N matrix
    EpetraMatrix(uint M, uint N);

    /// Copy constuctor
    explicit EpetraMatrix(const EpetraMatrix& A);

    /// Create matrix from given Epetra_FECrsMatrix pointer
    explicit EpetraMatrix(boost::shared_ptr<Epetra_FECrsMatrix> A);

    /// Create matrix from given Epetra_CrsGraph
    explicit EpetraMatrix(const Epetra_CrsGraph& graph);

    /// Destructor
    virtual ~EpetraMatrix();

    //--- Implementation of the GenericTensor interface ---

    /// Initialize zero tensor using sparsity pattern
    virtual void init(const GenericSparsityPattern& sparsity_pattern);

    /// Return copy of tensor
    virtual EpetraMatrix* copy() const;

    /// Return size of given dimension
    virtual uint size(uint dim) const;

    /// Set all entries to zero and keep any sparse structure
    virtual void zero();

    /// Finalize assembly of tensor
    virtual void apply();

    /// Display tensor
    virtual void disp(uint precision=2) const;

    //--- Implementation of the GenericMatrix interface ---

    /// Resize matrix to M x N
    virtual void resize(uint M, uint N);

    /// Get block of values
    virtual void get(double* block, uint m, const uint* rows, uint n, const uint* cols) const;

    /// Set block of values
    virtual void set(const double* block, uint m, const uint* rows, uint n, const uint* cols);

    /// Add block of values
    virtual void add(const double* block, uint m, const uint* rows, uint n, const uint* cols);

    /// Add multiple of given matrix (AXPY operation)
    virtual void axpy(double a, const GenericMatrix& A);

    /// Get non-zero values of given row
    virtual void getrow(uint row, Array<uint>& columns, Array<double>& values) const;

    /// Set values for given row
    virtual void setrow(uint row, const Array<uint>& columns, const Array<double>& values);

    /// Set given rows to zero
    virtual void zero(uint m, const uint* rows);

    /// Set given rows to identity matrix
    virtual void ident(uint m, const uint* rows);

    // Matrix-vector product, y = Ax
    virtual void mult(const GenericVector& x, GenericVector& y, bool transposed=false) const;

    /// Multiply matrix by given number
    virtual const EpetraMatrix& operator*= (double a);

    /// Divide matrix by given number
    virtual const EpetraMatrix& operator/= (double a);

    /// Assignment operator
    virtual const GenericMatrix& operator= (const GenericMatrix& x);
    //{ error("Not implemented."); return *this; }

    //--- Special functions ---

    /// Return linear algebra backend factory
    virtual LinearAlgebraFactory& factory() const;

    //--- Special Epetra functions ---

    /// Return Epetra_FECrsMatrix pointer
    boost::shared_ptr<Epetra_FECrsMatrix> mat() const;

    /// Assignment operator
    const EpetraMatrix& operator= (const EpetraMatrix& x);
    //{ error("Not implemented."); return *this; }

  private:

    // Epetra_FECrsMatrix pointer
    boost::shared_ptr<Epetra_FECrsMatrix> A;
    
  };

  LogStream& operator<< (LogStream& stream, const EpetraMatrix& A);

}

#endif //HAS_TRILINOS
#endif //__EPETRA_MATRIX_H
