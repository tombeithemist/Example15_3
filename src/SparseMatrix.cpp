/*
 * SparseMatrix.cpp
 *
 *  Created on: 22 Oct 2014
 *      Author: david
 */

#include <fstream>
#include <iostream>
#include <math.h>
#include <string>


#include "SparseMatrix.h"

using namespace std;

// Class definitions
class SparseMatrix {
public:
  SparseMatrix(int nrow, int nnzperrow) {
    // This constructor is called if we happen to know the number of rows
    // and an estimate of the number of nonzero entries per row.
    this->initialize(nrow, nnzperrow);
  }
  SparseMatrix() {
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
  ~SparseMatrix() {
    if (val_)
      delete[] val_;
    if (col_)
      delete[] col_;
    if (row_)
      delete[] row_;
    if (nnzs_)
      delete[] nnzs_;
  }
  void initialize(int nrow, int nnzperrow) {
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
  void finalize() {
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
  inline
  double& operator()(int m, int n) {
    // If the arrays are already full and inserting this entry would cause us to run off the end,
    // then we'll need to resize the arrays before inserting it
    if (nnzs_[m] >= N_nz_rowmax_) {
      this->reallocate();
    }
    // Search between row(m) and row(m+1) for col(k) = n (i.e. is the entry already in the matrix)
    int k = row_[m];
    bool foundEntry = false;
    while (k < (row_[m] + nnzs_[m]) && !foundEntry) {
      if (col_[k] == n) {
        foundEntry = true;
      }
      k++;
    }
    // If the entry is already in the matrix, then return a reference to it
    if (foundEntry) {
      return val_[k - 1];
    }
    // If the entry is not already in the matrix then we'll need to insert it
    else {
      N_nz_++;
      nnzs_[m]++;
      col_[k] = n;
      return val_[k];
    }
  }
  inline
  double& operator()(int k) {
    return val_[k];
  }
  void operator=(const SparseMatrix& A) {
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
  inline
  void multiply(double* u, double* v) {
    // Note: This function will perform a matrix vector multiplication with the input vector v, returning the output in u.
    for (int m = 0; m < N_row_; m++) {
      u[m] = 0.0;
      for (int k = row_[m]; k < row_[m + 1]; k++) {
        u[m] += val_[k] * v[col_[k]];
      }
    }
    return;
  }
  inline
  void multiply(double* u, double* v, bool* includerows, bool* includecols) {
    // Note: This function will perform a matrix vector multiplication on part of the matrix
    for (int m = 0; m < N_row_; m++) {
      u[m] = 0.0;
      if (includerows[m]) {
        for (int k = row_[m]; k < row_[m + 1]; k++) {

          if (includecols[col_[k]]) {
            u[m] += val_[k] * v[col_[k]];
          }
        }
      }
    }
    return;
  }
  inline
  void subtract(double u, SparseMatrix& A) {
    for (int k = 0; k < N_nz_; k++) {
      val_[k] -= (u * A.val_[k]);
    }
    return;
  }
  inline
  int getNnz() {
    return N_nz_;
  }
  inline
  int getNrow() {
    return N_row_;
  }
  void print(const char* name) {
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
protected:
  void reallocate() {
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
private:
  double* val_;         // [N_nz]    Stores the nonzero elements of the matrix
  int* col_;    // [N_nz]    Stores the column indices of the elements in each row
  int* row_;            // [N_row+1] Stores the locations in val that start a row
  int* nnzs_;   // [N_row+1] Stores the number of nonzero entries per row during the assembly process
  int N_row_;           // The number of rows in the matric
  int N_nz_;        // The number of non-zero entries currently stored in the matrix
  int N_nz_rowmax_; // The maximum number of non-zero entries per row. This will be an estimate until the matrix is assembled
  int N_allocated_; // The number of non-zero entries currently allocated for in val_ and col_
};
