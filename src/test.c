#include "header.h"

SEXP test(){

  Rprintf("jeg *har* skrevet en C-funksjon\n");
  Rprintf("jeg allokerer en såkalt *array* av floats\n");
  
  SEXP maxvalSEXP = PROTECT(allocVector(REALSXP, 10000000));
  
  double * array = REAL(maxvalSEXP); //pointer to array

  
  UNPROTECT(1);
  
  
  return NULL;
}

