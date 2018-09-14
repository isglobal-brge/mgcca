#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

// INPUT:
// Q -> matriz cuadrada (nxn)
// L -> vector (n)
// lambda -> escalar
// N

SEXP solve_penal(SEXP Q, SEXP L, SEXP LAMBDA, SEXP N){
  // Definicion de variables input y auxiliares
  int i, j, k, *n = INTEGER(N);
  double *q = REAL(Q), res_i, *res, *lambda = REAL(LAMBDA), *l = REAL(L);
  SEXP RES;
  
  // Definicion de output
  PROTECT(RES = allocMatrix(REALSXP, *n, *n));
  res = REAL(RES);
  
  // Calculos
  for(i = 0; i < *n; i++){
    for(j = i; j < *n; j++){
      res_i = 0;
      for(k = 0; k < *n; k++)
	//	Rprintf("%d %d, ", i
	res_i += (q[i + *n * k]*q[j + *n * k])/(l[k] + *lambda);
      res[j + i * *n] = res_i;
      res[i + j * *n] = res_i;
    }
  }

  // Unprotect y devolucion
  UNPROTECT(1);
  return(RES);
}
