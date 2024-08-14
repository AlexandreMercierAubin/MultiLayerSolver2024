#ifndef SOLVE_3BY3
#define SOLVE_3BY3

static void invA(double* A, double* Ainv)
{
    double det = A[0] * (A[4] * A[8] - A[5] * A[7]) -
         A[3] * (A[1] * A[8] - A[7] * A[2]) +
         A[6] * (A[1] * A[5] - A[4] * A[2]);
    
    double invdet = 1 / det;
    Ainv[0] = (A[4] * A[8] - A[5] * A[7]) * invdet;
    Ainv[1] = (A[6] * A[5] - A[3] * A[8]) * invdet;
    Ainv[2] = (A[3] * A[7] - A[6] * A[4]) * invdet;
    Ainv[3] = (A[7] * A[2] - A[1] * A[8]) * invdet;
    Ainv[4] = (A[0] * A[8] - A[6] * A[2]) * invdet;
    Ainv[5] = (A[1] * A[6] - A[0] * A[7]) * invdet;
    Ainv[6] = (A[1] * A[5] - A[2] * A[4]) * invdet;
    Ainv[7] = (A[2] * A[3] - A[0] * A[5]) * invdet;
    Ainv[8] = (A[0] * A[4] - A[1] * A[3]) * invdet;
    return;
} 

static void multMatVec(double* A, double* b, double* x)
{
    x[0] = A[0]*b[0]+A[3]*b[1]+A[6]*b[2];
    x[1] = A[1]*b[0]+A[4]*b[1]+A[7]*b[2];
    x[2] = A[2]*b[0]+A[5]*b[1]+A[8]*b[2];
    return;
} 

static void solve3by3(double* A, double* b, double* x)
{
  double* Ainv = new double[9];
  invA(A,Ainv);
  multMatVec(Ainv, b, x);
  delete[] Ainv;
  return;
} 

#endif