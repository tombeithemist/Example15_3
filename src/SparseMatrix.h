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
protected:
  void reallocate();

};

#endif /* SPARSEMATRIX_H_ */

