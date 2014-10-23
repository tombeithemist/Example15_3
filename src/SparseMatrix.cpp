/*
 * SparseMatrix.cpp
 *
 *  Created on: 22 Oct 2014
 *      Author: david
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>

#include "SparseMatrix.h"

using namespace std;

SparseMatrix::SparseMatrix(int nrow, int nnzperrow) {
  // This constructor is called if we happen to know the number of rows
  // and an estimate of the number of nonzero entries per row.
  this->initialize(nrow, nnzperrow);
}

SparseMatrix::SparseMatrix() {
  // This constructor is called if we have no useful information
  N_row_ = 0;
  N_nz_ = 0;
  N_nz_rowmax_ = 0;
  N_allocated_ = 0;
  val_ = NULL;
  col_ = NULL;
  row_ = NULL;
  nnzs_ = NULL;
}
SparseMatrix::~SparseMatrix() {
  if (val_)
    delete[] val_;
  if (col_)
    delete[] col_;
  if (row_)
    delete[] row_;
  if (nnzs_)
    delete[] nnzs_;
}

void SparseMatrix::initialize(int nrow, int nnzperrow) {
  N_row_ = nrow;
  N_nz_ = 0;
  N_nz_rowmax_ = nnzperrow;
  N_allocated_ = N_row_ * N_nz_rowmax_;
  val_ = new double[N_allocated_];
  col_ = new int[N_allocated_];
  row_ = new int[N_row_ + 1];
  nnzs_ = new int[N_row_];

  memset(val_, 0, N_allocated_ * sizeof(double));
  memset(col_, -1, N_allocated_ * sizeof(int));
  memset(row_, 0, (N_row_ + 1) * sizeof(int));
  memset(nnzs_, 0, (N_row_ + 1) * sizeof(int));

  for (int k = 0, kk = 0; k < N_row_; k++, kk += N_nz_rowmax_) {
    row_[k] = kk;
  }
  return;
}
void SparseMatrix::finalize() {
  int minCol = 0;
  int insertPos = 0;
  int index = 0;

  // Now that the matrix is assembled we can set N_nz_rowmax_ explicitly by
  // taking the largest value in the nnzs_ array
  N_nz_rowmax_ = 0;
  for (int m = 0; m < N_row_; m++) {
    N_nz_rowmax_ = max(N_nz_rowmax_, nnzs_[m]);
  }

  double* tempVal = new double[N_nz_];
  int* tempCol = new int[N_nz_];
  int* tempRow = new int[N_row_ + 1];
  bool* isSorted = new bool[N_allocated_]; // This array will help us sort the column indices

  memset(tempVal, 0, N_nz_ * sizeof(double));
  memset(tempCol, 0, N_nz_ * sizeof(int));
  memset(tempRow, 0, (N_row_ + 1) * sizeof(int));
  memset(isSorted, 0, N_allocated_ * sizeof(bool));

  for (int m = 0; m < N_row_; m++) {
    for (int k = row_[m]; k < (row_[m] + nnzs_[m]); k++) {
      minCol = N_row_ + 1;
      for (int kk = row_[m]; kk < (row_[m] + nnzs_[m]); kk++) {
        if (!isSorted[kk] && col_[kk] < minCol) {
          index = kk;
          minCol = col_[index];
        }
      }
      tempVal[insertPos] = val_[index];
      tempCol[insertPos] = col_[index];
      isSorted[index] = true;
      insertPos++;
    }
    tempRow[m + 1] = tempRow[m] + nnzs_[m];
  }

  delete[] val_;
  delete[] col_;
  delete[] row_;
  delete[] nnzs_;
  delete[] isSorted;

  val_ = tempVal;
  col_ = tempCol;
  row_ = tempRow;
  nnzs_ = NULL;
  N_allocated_ = N_nz_;

  return;
}

void SparseMatrix::operator=(const SparseMatrix& A) {
  if (val_)
    delete[] val_;
  if (col_)
    delete[] col_;
  if (row_)
    delete[] row_;
  if (nnzs_)
    delete[] nnzs_;

  N_row_ = A.N_row_;
  N_nz_ = A.N_nz_;
  N_nz_rowmax_ = A.N_nz_rowmax_;
  N_allocated_ = A.N_allocated_;
  val_ = new double[N_allocated_];
  col_ = new int[N_allocated_];
  row_ = new int[N_row_ + 1];

  memcpy(val_, A.val_, N_nz_ * sizeof(double));
  memcpy(col_, A.col_, N_nz_ * sizeof(int));
  memcpy(row_, A.row_, (N_row_ + 1) * sizeof(int));
}


void SparseMatrix::print(const char* name) {
  fstream matrix;
  cout << "Matrix " << name << " has " << N_row_ << " rows with " << N_nz_
      << " non-zero entries - " << N_allocated_ << " allocated." << flush;
  matrix.open(name, ios::out);
  matrix << "Mat = [" << endl;
  for (int m = 0; m < N_row_; m++) {
    for (int n = row_[m]; n < row_[m + 1]; n++) {
      matrix << m + 1 << "\t" << col_[n] + 1 << "\t" << val_[n] << endl;
    }
  }
  matrix << "];" << endl;
  matrix.close();
  cout << " Done." << flush << endl;
  return;
}

void SparseMatrix::reallocate() {
  // Double the memory allocation size
  N_nz_rowmax_ *= 2;

  N_allocated_ = N_nz_rowmax_ * N_row_;

  // Create some temporary arrays of the new size
  double* tempVal = new double[N_allocated_];
  int* tempCol = new int[N_allocated_];

  memset(tempVal, 0, N_allocated_ * sizeof(double));
  memset(tempCol, 0, N_allocated_ * sizeof(int));

  for (int m = 0, mm = 0; m < N_row_; m++, mm += N_nz_rowmax_) {
    memcpy(&tempVal[mm], &val_[row_[m]], nnzs_[m] * sizeof(double));
    memcpy(&tempCol[mm], &col_[row_[m]], nnzs_[m] * sizeof(int));
    row_[m] = mm;
  }

  // Delete the memory used by the old arrays
  delete[] val_;
  delete[] col_;

  // Assign the addresses of the new arrays
  val_ = tempVal;
  col_ = tempCol;

  return;
}

