/////////////////////////////////////////////////////////////////////////////
//
// Applied Numerical Methods
//
// Example:		15.3
//
// Problem:		dphi/dt + v . Grad phi = mu Grad^2 phi + psi
//
// Method:		Finite Element Method with linear 2D triangular elements
//				and Implicit Euler Method and Conjugate Gradient Method
//
// Compilation:	g++ Example15_3.cpp -o Example15_3
//
// Execution:	./Example15_3 Example15_3.grid
//
/////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <math.h>

#include "SparseMatrix.h"
#include "Boundary.h"


using namespace std;


// Global variables
const double t_min = 0.00;
const double t_max = 2.00;
const double Delta_t = 0.01;
const double v[2] = { 0.50, 0.50 };
const double mu = 0.01;
const double psi = 0.20;
const int N_t = static_cast<int>((t_max - t_min) / Delta_t + 1);

// Function declarations
int read(char* filename, double**& Points, int**& Faces, int**& Elements,
    Boundary*& Boundaries, int& N_p, int& N_f, int& N_e, int& N_b);
void write(fstream& file, double* phi, int N_p);
void assemble(SparseMatrix& M, SparseMatrix& K, double* s, double* phi,
    bool* Free, bool* Fixed, double** Points, int** Faces, int** Elements,
    Boundary* Boundaries, int N_p, int N_f, int N_e, int N_b);
void solve(SparseMatrix& A, double* phi, double* b, bool* Free, bool* Fixed);

int main(int argc, char** argv) {
  // Simulation parameters
  double** Points = NULL;
  int** Faces = NULL;
  int** Elements = NULL;
  Boundary* Boundaries = NULL;
  int N_p = 0;
  int N_f = 0;
  int N_e = 0;
  int N_b = 0;
  double t = 0;
  fstream file;

  if (argc < 2) {
    cerr << "No grid file specified" << endl;
    return 1;
  } else {
    int rc = read(argv[1], Points, Faces, Elements, Boundaries, N_p, N_f, N_e, N_b);
    if(rc != 0) {
      cerr << "Exiting as a result of error. rc: " << rc << endl;
      return(rc);
    }

  }
  // Allocate arrays
  double* phi = new double[N_p];
  double* s = new double[N_p];
  double* b = new double[N_p];
  bool* Free = new bool[N_p];
  bool* Fixed = new bool[N_p];
  double* AphiFixed = new double[N_p];
  SparseMatrix M;
  SparseMatrix K;
  SparseMatrix A;

  // Set initial condition
  t = t_min;
  for (int m = 0; m < N_p; m++) {
    phi[m] = exp(
        -50 * (pow(Points[m][0] - 0.3, 2) + pow(Points[m][1] - 0.3, 2))) + 1.0;
  }

  assemble(M, K, s, phi, Free, Fixed, Points, Faces, Elements, Boundaries, N_p,
      N_f, N_e, N_b);

  A = M;
  A.subtract(Delta_t, K); // At this point we have A = M-Delta_t*K

  // Compute the column vector to subtract from the right hand side to take account of fixed nodes
  A.multiply(AphiFixed, phi, Free, Fixed);

  file.open("Example15_3.data", ios::out);
  write(file, phi, N_p);

  // Time marching loop
  for (int l = 0; l < N_t - 1; l++) {
    t += Delta_t;
    cout << "t = " << t;

    // Assemble b
    M.multiply(b, phi);
    for (int m = 0; m < N_p; m++) {
      b[m] += Delta_t * s[m] - AphiFixed[m];
    }

    // Solve the linear system
    solve(A, phi, b, Free, Fixed);

    // Write the solution
    write(file, phi, N_p);
  }

  file.close();

  // Deallocate arrays
  for (int boundary = 0; boundary < N_b; boundary++) {
    delete[] Boundaries[boundary].indices_;
  }
  delete[] Points[0];
  delete[] Points;
  delete[] Faces[0];
  delete[] Faces;
  delete[] Elements[0];
  delete[] Elements;
  delete[] Boundaries;
  delete[] phi;
  delete[] s;
  delete[] b;
  delete[] Free;
  delete[] Fixed;
  delete[] AphiFixed;

  return 0;
}

int read(char* filename, double**& Points, int**& Faces, int**& Elements,
    Boundary*& Boundaries, int& N_p, int& N_f, int& N_e, int& N_b) {
  fstream file;
  string temp;

  cout << "Reading " << filename << "... " << flush;

  file.open(filename);
  if (!file.is_open()) {
    cerr << "Error opening file: " << filename << endl;
    return(1);
  }

  file >> temp >> N_p;
  file >> temp >> N_f;
  file >> temp >> N_e;
  file >> temp >> N_b;

  Points = new double*[N_p];
  Faces = new int*[N_f];
  Elements = new int*[N_e];
  Boundaries = new Boundary[N_b];
  Points[0] = new double[N_p * 2];
  Faces[0] = new int[N_f * 2];
  Elements[0] = new int[N_e * 3];
  for (int p = 1, pp = 2; p < N_p; p++, pp += 2) {
    Points[p] = &Points[0][pp];
  }
  for (int f = 1, ff = 2; f < N_f; f++, ff += 2) {
    Faces[f] = &Faces[0][ff];
  }
  for (int e = 1, ee = 3; e < N_e; e++, ee += 3) {
    Elements[e] = &Elements[0][ee];
  }

  file >> temp;
  for (int p = 0; p < N_p; p++) {
    file >> Points[p][0] >> Points[p][1];
  }

  file >> temp;
  for (int f = 0; f < N_f; f++) {
    file >> Faces[f][0] >> Faces[f][1];
  }

  file >> temp;
  for (int e = 0; e < N_e; e++) {
    file >> Elements[e][0] >> Elements[e][1] >> Elements[e][2];
  }

  file >> temp;
  for (int b = 0; b < N_b; b++) {
    file >> Boundaries[b].name_ >> Boundaries[b].type_ >> Boundaries[b].N_;
    Boundaries[b].indices_ = new int[Boundaries[b].N_];
    for (int n = 0; n < Boundaries[b].N_; n++) {
      file >> Boundaries[b].indices_[n];
    }
    file >> Boundaries[b].value_;
  }

  file.close();

  cout << "Done.\n" << flush;

  return(0);
}

void write(fstream& file, double* phi, int N_p) {
  for (int m = 0; m < N_p; m++) {
    file << phi[m] << "\t";
  }
  file << "\n";
  return;
}

void assemble(SparseMatrix& M, SparseMatrix& K, double* s, double* phi,
    bool* Free, bool* Fixed, double** Points, int** Faces, int** Elements,
    Boundary* Boundaries, int N_p, int N_f, int N_e, int N_b) {
  cout << "Assembling system... " << flush;

  double x[3];
  double y[3];
  double gradEta[2][3];
  double gradEta_p[2] = { 0.0, 0.0 };
  double gradEta_q[2] = { 0.0, 0.0 };
  double M_e[3][3] = { { 2.0, 1.0, 1.0 }, { 1.0, 2.0, 1.0 }, { 1.0, 1.0, 2.0 } };
  double s_e[3] = { 1.0, 1.0, 1.0 };
  int Nodes[3] = { 0, 0, 0 };
  double* Omega = new double[N_e];
  double* Gamma = new double[N_f];
  int m;
  int n;

  // Assign all the indices to be free initially
  for (int p = 0; p < N_p; p++) {
    Free[p] = 1;
    Fixed[p] = 0;
    s[p] = 0.0;
  }

  // Calculate face lengths
  for (int f = 0; f < N_f; f++) {
    for (int p = 0; p < 2; p++) {
      x[p] = Points[Faces[f][p]][0];
      y[p] = Points[Faces[f][p]][1];
    }
    Gamma[f] = sqrt(
        (x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]));
  }

  // Calculate element areas
  for (int e = 0; e < N_e; e++) {
    for (int p = 0; p < 3; p++) {
      x[p] = Points[Elements[e][p]][0];
      y[p] = Points[Elements[e][p]][1];
    }
    Omega[e] = ((x[0] * y[1] - x[1] * y[0]) + (x[2] * y[0]) - (x[0] * y[2])
        + (x[1] * y[2] - x[2] * y[1])) / 2;
  }

  // Assemble M, K, and s
  M.initialize(N_p, 10);
  K.initialize(N_p, 10);

  for (int e = 0; e < N_e; e++) {
    for (int p = 0; p < 3; p++) {
      Nodes[p] = Elements[e][p];
      x[p] = Points[Nodes[p]][0];
      y[p] = Points[Nodes[p]][1];
    }

    gradEta[0][0] = (y[1] - y[2]) / (2 * Omega[e]);
    gradEta[0][1] = (y[2] - y[0]) / (2 * Omega[e]);
    gradEta[0][2] = (y[0] - y[1]) / (2 * Omega[e]);
    gradEta[1][0] = (x[2] - x[1]) / (2 * Omega[e]);
    gradEta[1][1] = (x[0] - x[2]) / (2 * Omega[e]);
    gradEta[1][2] = (x[1] - x[0]) / (2 * Omega[e]);

    // Outer loop over each node
    for (int p = 0; p < 3; p++) {
      m = Nodes[p];
      gradEta_p[0] = gradEta[0][p];
      gradEta_p[1] = gradEta[1][p];

      // Inner loop over each node
      for (int q = 0; q < 3; q++) {
        n = Nodes[q];
        gradEta_q[0] = gradEta[0][q];
        gradEta_q[1] = gradEta[1][q];

        M(m, n) += M_e[p][q] * Omega[e] / 12;

        // There is a missing right parenthesis in this line???
        // K(m,n) -= ((v[0]*gradEta_q[0]+v[1]*gradEta_q[1])/3*Omega[e] + mu*(gradEta_p[0]*gradEta_q[0]+gradEta_p[1]*gradEta_q[1])*Omega[e];
        K(m, n) -= (v[0] * gradEta_q[0] + v[1] * gradEta_q[1]) / 3 * Omega[e]
            + mu * (gradEta_p[0] * gradEta_q[0] + gradEta_p[1] * gradEta_q[1])
                * Omega[e];
      }
      s[m] += s_e[p] * psi * Omega[e] / 3;
    }
  }

  // Apply boundary conditions
  for (int b = 0; b < N_b; b++) {
    if (Boundaries[b].type_ == "neumann") {
      for (int f = 0; f < Boundaries[b].N_; f++) {
        for (int p = 0; p < 2; p++) {
          Nodes[p] = Faces[Boundaries[b].indices_[f]][p];
          m = Nodes[p];
          s[m] += mu * Boundaries[b].value_ * Gamma[Boundaries[b].indices_[f]]
              / 2;
        }
      }
    } else if (Boundaries[b].type_ == "dirichlet") {
      for (int p = 0; p < Boundaries[b].N_; p++) {
        m = Boundaries[b].indices_[p];
        phi[m] = Boundaries[b].value_;
        Free[m] = false;
        Fixed[m] = true;
      }
    }
  }

  K.finalize();
  M.finalize();

  delete[] Gamma;
  delete[] Omega;

  cout << "Done.\n" << flush;

  return;
}

void solve(SparseMatrix& A, double* phi, double* b, bool* Free, bool* Fixed) {
  int N_row = A.getNrow();
  double* r_old = new double[N_row];
  double* r = new double[N_row];
  double* d = new double[N_row];
  double* Ad = new double[N_row];
  double* Aphi = new double[N_row];
  double alpha = 0.0;
  double beta = 0.0;
  double r_norm = 0.0;
  double tolerance = 1e-8;
  double N_k = 1e+3;
  double r_oldTr_old = 0.0;
  double rTr = 0.0;
  double dTAd = 0.0;
  int k = 0;
  int m = 0;
  int n = 0;

  memset(r_old, 0, N_row * sizeof(double));
  memset(r, 0, N_row * sizeof(double));
  memset(d, 0, N_row * sizeof(double));
  memset(Ad, 0, N_row * sizeof(double));

  // Compute the initial residual
  A.multiply(Aphi, phi, Free, Free);
  for (m = 0; m < N_row; m++) {
    if (Free[m]) {
      r_old[m] = b[m] - Aphi[m];
      d[m] = r_old[m];
      r_oldTr_old += r_old[m] * r_old[m];
    }
  }
  r_norm = sqrt(r_oldTr_old);

  // Conjugate Gradient iterative loop
  while (r_norm > tolerance && k < N_k) {
    dTAd = 0.0;
    A.multiply(Ad, d, Free, Free);
    for (m = 0; m < N_row; m++) {
      if (Free[m]) {
        dTAd += d[m] * Ad[m];
      }
    }
    alpha = r_oldTr_old / dTAd;
    for (m = 0; m < N_row; m++) {
      if (Free[m]) {
        phi[m] += alpha * d[m];
      }
    }
    for (m = 0; m < N_row; m++) {
      if (Free[m]) {
        r[m] = r_old[m] - alpha * Ad[m];
      }
    }
    rTr = 0.0;
    for (m = 0; m < N_row; m++) {
      if (Free[m]) {
        rTr += r[m] * r[m];
      }
    }
    beta = rTr / r_oldTr_old;
    for (m = 0; m < N_row; m++) {
      if (Free[m]) {
        d[m] = r[m] + beta * d[m];
      }
    }
    for (m = 0; m < N_row; m++) {
      if (Free[m]) {
        r_old[m] = r[m];
      }
    }
    r_oldTr_old = rTr;
    r_norm = sqrt(rTr);
    k++;
  }

  cout << ", k = " << k << ", r_norm = " << r_norm << endl;

  delete[] r_old;
  delete[] r;
  delete[] d;
  delete[] Ad;
  delete[] Aphi;

  return;
}

