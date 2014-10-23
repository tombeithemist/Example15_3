/*
 * SparseMatrix.h
 *
 *  Created on: 22 Oct 2014
 *      Author: david
 */


#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

class SparseMatrix {
  public:
    SparseMatrix(int nrow, int nnzperrow);
    SparseMatrix();
    virtual ~SparseMatrix();
    void initialize(int nrow, int nnzperrow);
    void finalize();
    double& operator()(int m, int n);
    double& operator()(int k);
    void operator=(const SparseMatrix& A);
    void multiply(double* u, double* v);
    void multiply(double* u, double* v, bool* includerows, bool* includecols);
    void subtract(double u, SparseMatrix& A);
    int getNnz();
    int getNrow();
    void print(const char* name);
  protected:
    void reallocate();
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

// ----------------------------------------------------------------------------------------
// Definition of inline methods needs to be in header file so they are only included once.
// ----------------------------------------------------------------------------------------

inline
double& SparseMatrix::operator()(int m, int n) {
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
double& SparseMatrix::operator()(int k) {
  return val_[k];
}

inline
void SparseMatrix::multiply(double* u, double* v) {
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
void SparseMatrix::multiply(double* u, double* v, bool* includerows,
    bool* includecols) {
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
void SparseMatrix::subtract(double u, SparseMatrix& A) {
  for (int k = 0; k < N_nz_; k++) {
    val_[k] -= (u * A.val_[k]);
  }
  return;
}

inline
int SparseMatrix::getNnz() {
  return N_nz_;
}

inline
int SparseMatrix::getNrow() {
  return N_row_;
}

#endif /* SPARSEMATRIX_H_ */

