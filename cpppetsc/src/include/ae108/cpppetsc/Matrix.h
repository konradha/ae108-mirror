// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "ae108/cpppetsc/InsertionProxy.h"
#include "ae108/cpppetsc/InvalidParametersException.h"
#include "ae108/cpppetsc/Matrix_fwd.h"
#include "ae108/cpppetsc/Mesh_fwd.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include <cstring>
#include <initializer_list>
#include <memory>
#include <numeric>
#include <petscdm.h>
#include <petscmat.h>
#include <petscmath.h>
#include <petscsys.h>
#include <utility>
#include <vector>

namespace ae108 {
namespace cpppetsc {

template <class Policy> class Matrix {
public:
  using size_type = PetscInt;
  using value_type = PetscScalar;
  using real_type = PetscReal;

  /**
   * @brief Allocates a matrix with global_rows rows and global_columns columns.
   */
  explicit Matrix(size_type global_rows, size_type global_columns);

  /**
   * @brief Creates a Matrix from a list of rows.
   *
   * @throw InvalidParametersException if the number of columns is inconsistent.
   */
  static Matrix
  fromList(const std::initializer_list<std::initializer_list<value_type>> list);

  /**
   * @brief Creates a matrix with a structure defined by the given mesh.
   */
  static Matrix fromMesh(const Mesh<Policy> &mesh);

  /**
   * @brief Creates a matrix with the same layout and nonzero pattern as the
   * provided matrix.
   */
  static Matrix fromLayoutOf(const Matrix &matrix);

  /**
   * @brief Creates a matrix A * B * C.
   */
  static Matrix fromProduct(const Matrix &A, const Matrix &B, const Matrix &C);

  /**
   * @brief Creates a matrix A * B.
   */
  static Matrix fromProduct(const Matrix &A, const Matrix &B);

  /**
   * @brief creates a matrix A^t * B.
   */

  static Matrix fromAtB(const Matrix &A, const Matrix &B);

  /**
   * @brief Creates a matrix P^t * A * P.
   */
  static Matrix fromPtAP(const Matrix &P, const Matrix &A);

  /**
   * @brief Creates a matrix P * A * P^t.
   */
  static Matrix fromPAPt(const Matrix &P, const Matrix &A);

  /**
   * @brief Creates a matrix from the inverted block diagonal.
   */
  static Matrix fromInvertedBlockDiagonalOf(const Matrix &matrix);

  /**
   * @brief Returns a deep copy of the matrix.
   */
  static Matrix clone(const Matrix &matrix);

  /**
   * @brief Wraps the given UniqueEntity in a Matrix.
   */
  explicit Matrix(UniqueEntity<Mat> mat);

  /**
   * @brief Returns a pair of (number of rows, number of columns).
   */
  std::pair<size_type, size_type> size() const;

  /**
   * @brief Returns the block size of the matrix.
   */
  size_type blockSize() const;

  class AssemblyView;
  /**
   * @brief Returns a view that permits assembly of the matrix (e.g. adding
   * entries to the current entries). If the number of nonzero entries per row
   * is known in advance, then it is more efficient to use the method
   * preallocatedAssemblyView(...).
   *
   * @remark Final assembly is automatically performed when the view and its
   * inserters have been destroyed.
   */
  AssemblyView assemblyView();

  /**
   * @brief Returns a view that permits assembly of the matrix (e.g. adding
   * entries to the current entries).
   *
   * @param nonZeroesPerRow The maximum number of nonzero values per row to
   * expect. If the matrix is of type MPIAIJ, then this parameter defines the
   * number of nonzero values in the diagonal block. No space for
   * values in the off-diagonal block is allocated in this case.
   *
   * @throw InvalidParametersException if matrix is not of type SEQAIJ or
   * MPIAIJ.
   *
   * @remark Final assembly is automatically performed when the view and its
   * inserters have been destroyed.
   */
  AssemblyView preallocatedAssemblyView(const size_type nonZeroesPerRow);

  /**
   * @brief Computes the Frobenius norm of the matrix.
   */
  real_type norm() const;

  /**
   * @brief Print the matrix to world stdout.
   */
  void print() const;

  /**
   * @brief Returns the local row indices in the format:
   * (first row index, last row index + 1).
   */
  std::pair<size_type, size_type> localRowRange() const;

  /**
   * @brief Returns the entry at the provided coordinates.
   *
   * @remark Does not check bounds.
   */
  value_type operator()(size_type row, size_type column) const;

  /**
   * @brief Adds alpha * X to the matrix.
   *
   * @pre The matrix X has the same nonzero pattern.
   */
  void addAlphaX(const value_type alpha, const Matrix &X);

  /**
   * @brief Replaces row indices by the corresponding rows of the identity
   * matrix. Does not change the nonzero pattern.
   *
   * @param rows Contains the global row indices to replace.
   *
   * @remark Must be called on every process in the parallel setting.
   */
  void replaceRowsByEye(const std::vector<size_type> &rows);

  /**
   * @brief Replaces all entries by zero, but keeps the structure for sparse
   * matrices.
   */
  void setZero();

  /**
   * @brief Scale the matrix by factor.
   */
  void scale(const value_type factor);

  /**
   * @brief Perform final assembly of the matrix.
   */
  void finalize();

  /**
   * @brief Returns the internal matrix.
   *
   * @remark For internal use only.
   */
  Mat data() const;

private:
  static Mat createMat();

  UniqueEntity<Mat> _mat;
};

template <class Policy> class Matrix<Policy>::AssemblyView {
public:
  /**
   * @param matrix nonzero pointer to a Matrix
   */
  explicit AssemblyView(std::shared_ptr<Matrix> matrix);

  /**
   * @brief Provides element() and elements() access functions to set/add
   * single/multiple elements. These calls can be chained (matrix is a Matrix
   * object):
   *
   * \code{.cpp}matrix.replace().element(0, 0, 7.).element(1, 1, 3.);\endcode
   *
   * It is also possible to set/add entries by accessing the index first:
   *
   * \code{.cpp}matrix.replace()(0, 0) = 7.;\endcode
   *
   * This way, the Inserter class can be used as a "matrix" for assembly.
   */
  template <InsertMode Mode> class Inserter;

  /**
   * @brief Replace elements with different elements. See the documentation of
   * the Inserter class for more information on how to use the return value.
   *
   * @remark Automatically broadcasts changes. Do not use add and replace
   * Inserters at the same time.
   */
  Inserter<INSERT_VALUES> replace() const;

  /**
   * @brief Add elements to current elements. See the documentation of the
   * Inserter class for more information on how to use the return value.
   *
   * @remark Automatically broadcasts changes. Do not use add and replace
   * Inserters at the same time.
   */
  Inserter<ADD_VALUES> add() const;

private:
  std::shared_ptr<Matrix> _matrix;
};

extern template class Matrix<SequentialComputePolicy>;
extern template class Matrix<ParallelComputePolicy>;
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

#include "ae108/cpppetsc/Mesh.h"

namespace ae108 {
namespace cpppetsc {

template <class Policy> Mat Matrix<Policy>::createMat() {
  Mat mat;
  Policy::handleError(MatCreate(Policy::communicator(), &mat));
  return mat;
}

template <class Policy>
Matrix<Policy> Matrix<Policy>::fromMesh(const Mesh<Policy> &mesh) {
  auto mat = Mat();
  Policy::handleError(DMCreateMatrix(mesh.data(), &mat));
  auto matrix = Matrix(makeUniqueEntity<Policy>(mat));

  return matrix;
}

template <class Policy>
Matrix<Policy> Matrix<Policy>::fromLayoutOf(const Matrix &matrix) {
  auto mat = Mat();
  Policy::handleError(
      MatDuplicate(matrix.data(), MAT_DO_NOT_COPY_VALUES, &mat));
  auto result = Matrix(makeUniqueEntity<Policy>(mat));

  return result;
}

template <class Policy>
Matrix<Policy> Matrix<Policy>::fromProduct(const Matrix &A, const Matrix &B,
                                           const Matrix &C) {
  auto mat = Mat();
  Policy::handleError(MatMatMatMult(A.data(), B.data(), C.data(),
                                    MAT_INITIAL_MATRIX, PETSC_DEFAULT, &mat));
  auto result = Matrix(makeUniqueEntity<Policy>(mat));

  return result;
}

template <class Policy>
Matrix<Policy> Matrix<Policy>::fromProduct(const Matrix &A, const Matrix &B) {
  auto mat = Mat();
  Policy::handleError(
      MatMatMult(A.data(), B.data(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &mat));
  auto result = Matrix(makeUniqueEntity<Policy>(mat));

  return result;
}

template <class Policy>
Matrix<Policy> Matrix<Policy>::fromAtB(const Matrix &A, const Matrix &B) {
  auto mat = Mat();
  Policy::handleError(MatTransposeMatMult(
      A.data(), B.data(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &mat));
  auto result = Matrix(makeUniqueEntity<Policy>(mat));

  return result;
}

template <class Policy>
Matrix<Policy> Matrix<Policy>::fromPtAP(const Matrix &P, const Matrix &A) {
  auto mat = Mat();
  Policy::handleError(
      MatPtAP(A.data(), P.data(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &mat));
  auto result = Matrix(makeUniqueEntity<Policy>(mat));

  return result;
}

template <class Policy>
Matrix<Policy> Matrix<Policy>::fromPAPt(const Matrix &P, const Matrix &A) {
  auto mat = Mat();
  Policy::handleError(
      MatRARt(A.data(), P.data(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &mat));
  auto result = Matrix(makeUniqueEntity<Policy>(mat));

  return result;
}

template <class Policy>
Matrix<Policy>
Matrix<Policy>::fromInvertedBlockDiagonalOf(const Matrix &matrix) {
  auto result = fromLayoutOf(matrix);

  Policy::handleError(MatInvertBlockDiagonalMat(matrix.data(), result.data()));
  return result;
}

template <class Policy>
Matrix<Policy> Matrix<Policy>::clone(const Matrix &matrix) {
  auto mat = Mat();
  Policy::handleError(MatDuplicate(matrix.data(), MAT_COPY_VALUES, &mat));
  auto result = Matrix(makeUniqueEntity<Policy>(mat));

  return result;
}

template <class Policy>
Matrix<Policy>::Matrix(UniqueEntity<Mat> mat) : _mat(std::move(mat)) {}

template <class Policy>
Matrix<Policy>::Matrix(const size_type global_rows,
                       const size_type global_columns)
    : _mat(makeUniqueEntity<Policy>(createMat())) {
  Policy::handleError(MatSetFromOptions(_mat.get()));
  Policy::handleError(MatSetSizes(_mat.get(), PETSC_DETERMINE, PETSC_DETERMINE,
                                  global_rows, global_columns));
  Policy::handleError(MatSetUp(_mat.get()));
  finalize();
}

template <class Policy>
Matrix<Policy> Matrix<Policy>::fromList(
    const std::initializer_list<std::initializer_list<value_type>> list) {
  auto result = Matrix(list.size(),
                       list.begin() == list.end() ? 0 : list.begin()->size());
  std::vector<size_type> column_indices(
      list.begin() == list.end() ? 0 : list.begin()->size());
  std::iota(column_indices.begin(), column_indices.end(), size_type{0});

  size_type row_index = 0;
  auto replacer = result.assemblyView().replace();
  for (const auto &rowlist : list) {
    if (rowlist.size() != column_indices.size())
      throw InvalidParametersException{};
    replacer.elements({row_index}, column_indices, rowlist);
    ++row_index;
  }
  return result;
}

template <class Policy>
typename Matrix<Policy>::AssemblyView Matrix<Policy>::assemblyView() {
  return AssemblyView(std::shared_ptr<Matrix>(
      this, [](Matrix *const matrix) { matrix->finalize(); }));
}

template <class Policy>
typename Matrix<Policy>::AssemblyView
Matrix<Policy>::preallocatedAssemblyView(const size_type nonZeroesPerRow) {
  auto type = MatType();
  Policy::handleError(MatGetType(_mat.get(), &type));

  if (std::strcmp(type, MATSEQAIJ) == 0) {
    MatSeqAIJSetPreallocation(_mat.get(), nonZeroesPerRow, nullptr);
  } else if (std::strcmp(type, MATMPIAIJ) == 0) {
    MatMPIAIJSetPreallocation(_mat.get(), nonZeroesPerRow, nullptr, 0, nullptr);
  } else {
    throw InvalidParametersException{};
  }

  return assemblyView();
}

template <class Policy>
Matrix<Policy>::AssemblyView::AssemblyView(std::shared_ptr<Matrix> matrix)
    : _matrix(matrix) {}

template <class Policy>
template <InsertMode Mode>
class Matrix<Policy>::AssemblyView::Inserter {

public:
  using value_type = Matrix::value_type;
  using size_type = Matrix::size_type;

  /**
   * @param matrix nonzero pointer to a Matrix
   */
  Inserter(std::shared_ptr<Matrix> matrix) : _matrix(matrix) {}

  size_type rows() const { return _matrix->size().first; }

  size_type cols() const { return _matrix->size().second; }

  const Inserter &elements(const std::vector<size_type> &row_indices,
                           const std::vector<size_type> &column_indices,
                           const std::vector<value_type> &values) const {
    Policy::handleError(MatSetValues(
        _matrix->_mat.get(), row_indices.size(), row_indices.data(),
        column_indices.size(), column_indices.data(), values.data(), Mode));
    return *this;
  }

  InsertionProxy<Mode> operator()(const size_type row_index,
                                  const size_type column_index) const {
    return InsertionProxy<Mode>{
        [this, row_index, column_index](const value_type value) {
          this->element(row_index, column_index, value);
        }};
  }

  const Inserter &element(const size_type row_index,
                          const size_type column_index,
                          const value_type value) const {
    Policy::handleError(
        MatSetValue(_matrix->_mat.get(), row_index, column_index, value, Mode));
    return *this;
  }

  void setZero() const {
    flush();
    _matrix->setZero();
  }

  ~Inserter() { flush(); }

private:
  void flush() const {
    Policy::handleError(
        MatAssemblyBegin(_matrix->_mat.get(), MAT_FLUSH_ASSEMBLY));
    Policy::handleError(
        MatAssemblyEnd(_matrix->_mat.get(), MAT_FLUSH_ASSEMBLY));
  }

  std::shared_ptr<Matrix> _matrix;
};

template <class Policy>
std::pair<typename Matrix<Policy>::size_type,
          typename Matrix<Policy>::size_type>
Matrix<Policy>::size() const {
  auto return_value = std::pair<size_type, size_type>();
  Policy::handleError(
      MatGetSize(_mat.get(), &return_value.first, &return_value.second));
  return return_value;
}

template <class Policy>
typename Matrix<Policy>::size_type Matrix<Policy>::blockSize() const {
  auto size = size_type{};
  Policy::handleError(MatGetBlockSize(_mat.get(), &size));
  return size;
}

template <class Policy>
std::pair<typename Matrix<Policy>::size_type,
          typename Matrix<Policy>::size_type>
Matrix<Policy>::localRowRange() const {
  auto return_value = std::pair<size_type, size_type>();
  Policy::handleError(MatGetOwnershipRange(_mat.get(), &return_value.first,
                                           &return_value.second));
  return return_value;
}

template <class Policy>
typename Matrix<Policy>::real_type Matrix<Policy>::norm() const {
  auto result = real_type{};
  Policy::handleError(MatNorm(_mat.get(), NORM_FROBENIUS, &result));
  return result;
}

template <class Policy> void Matrix<Policy>::print() const {
  Policy::handleError(MatView(_mat.get(), nullptr));
}

template <class Policy>
typename Matrix<Policy>::AssemblyView::template Inserter<INSERT_VALUES>
Matrix<Policy>::AssemblyView::replace() const {
  return Inserter<INSERT_VALUES>(_matrix);
}

template <class Policy>
typename Matrix<Policy>::AssemblyView::template Inserter<ADD_VALUES>
Matrix<Policy>::AssemblyView::add() const {
  return Inserter<ADD_VALUES>(_matrix);
}

template <class Policy>
typename Matrix<Policy>::value_type
Matrix<Policy>::operator()(const size_type row, const size_type column) const {
  const size_type row_indices[] = {row};
  const size_type column_indices[] = {column};
  value_type values[1] = {0.};
  Policy::handleError(
      MatGetValues(_mat.get(), 1, row_indices, 1, column_indices, values));
  return values[0];
}

template <class Policy>
void Matrix<Policy>::addAlphaX(const value_type alpha, const Matrix &X) {
  Policy::handleError(
      MatAXPY(_mat.get(), alpha, X.data(), SAME_NONZERO_PATTERN));
}

template <class Policy>
void Matrix<Policy>::replaceRowsByEye(const std::vector<size_type> &rows) {
  auto flag = PetscBool();
  Policy::handleError(
      MatGetOption(_mat.get(), MAT_KEEP_NONZERO_PATTERN, &flag));
  Policy::handleError(
      MatSetOption(_mat.get(), MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
  Policy::handleError(
      MatZeroRows(_mat.get(), rows.size(), rows.data(), 1.0, nullptr, nullptr));
  Policy::handleError(MatSetOption(_mat.get(), MAT_KEEP_NONZERO_PATTERN, flag));
}

template <class Policy> void Matrix<Policy>::setZero() {
  Policy::handleError(MatZeroEntries(_mat.get()));
}

template <class Policy> void Matrix<Policy>::scale(const value_type factor) {
  Policy::handleError(MatScale(_mat.get(), factor));
}

template <class Policy> void Matrix<Policy>::finalize() {
  Policy::handleError(MatAssemblyBegin(_mat.get(), MAT_FINAL_ASSEMBLY));
  Policy::handleError(MatAssemblyEnd(_mat.get(), MAT_FINAL_ASSEMBLY));
}

template <class Policy> Mat Matrix<Policy>::data() const { return _mat.get(); }
} // namespace cpppetsc
} // namespace ae108