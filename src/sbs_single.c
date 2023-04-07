#include "header.h"

SEXP cSBS_single(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI,SEXP rescale_variance_boolI,
    SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(thresholdI);
    
/*    PROTECT(asI);
    PROTECT(nu_asI);
*/    PROTECT(nI);
    PROTECT(pI);
/*    PROTECT(len_asI);
*/    PROTECT(debugI);
    PROTECT(rescale_variance_boolI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double threshold = *(REAL(thresholdI));
    
    int debug = *INTEGER(debugI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);

    //UNPROTECT(3); // unprotecting all except the arrays
    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("n = %d\n", n);
    }

    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vectorSEXP = PROTECT(allocVector(REALSXP, maxlen));

   
    double * vector = REAL(vectorSEXP);
    memset(vector, 0, sizeof(double)*maxlen);

    /*SEXP tmpvecSEXP = PROTECT(allocVector(REALSXP, len_as));
    double * tmpvec = REAL(tmpvecSEXP);
    memset(tmpvec, 0, sizeof(double)*len_as);*/

    if (rescale_variance_bool)
    {
        rescale_variance(X,NULL, n, p, vector);
    }

    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));

    for (int j = 1; j <=n; ++j)
    {
        for (int i = 0; i < p; ++i)
        {
            //#define cord_spec(r,c, D) ((r) + (D)*(c))
            cumsums[cord_spec(i,j,p)] = X[cord_spec(i,j-1, p)] +cumsums[cord_spec(i,j-1, p)];
        }
    }

    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p*(n)));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p*(n));






    //double prev_nu_a = as[0];
    //double a;
    //double nu_a;
    double tmp = 0;
    //double a_tmp=-100000;
    //int pos_a_tmp = 0;
    //double * val=0;
    //int n = e-s;

    //int detected = 0;
    double tmp2 =0;
    int s = -1;
    int e = n-1;
    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, 1));
    int * maxpos = INTEGER(maxposSEXP);
    *maxpos = 0;
    double maxval = NEGINF;

    CUSUM(cumsums, cusum, s, e, p);
    int end = e-s-1;
    for (int j = 0; j < end; ++j){
        /*memset(tmpvec, 0, sizeof(double)*len_as);*/
        tmp = 0;
        // first aggregate thresholded CUSUMs
        for (int i = 0; i < p; ++i){
            tmp2 = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
            
            if (tmp2> threshold)
            {
                tmp+= tmp2;

            
            }
        }
        if(tmp>maxval){
            maxval = tmp;
            *maxpos = j;
        }
    }

    SEXP maxvalSEXP = PROTECT(allocVector(REALSXP, 1));
    double * maxvall = REAL(maxvalSEXP);
    *maxvall = maxval;

    SEXP ret = PROTECT(allocVector(VECSXP, 2)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, maxposSEXP);
    SET_VECTOR_ELT(ret, 1, maxvalSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("pos"));
    SET_STRING_ELT(names, 1, mkChar("maxval"));


    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(13);
    return ret;
}





SEXP cSBS_single_calibrate(SEXP nI, SEXP pI, SEXP NI, SEXP tolnI, SEXP rescale_variance_boolI,SEXP debugI){

    PROTECT(nI);
    PROTECT(pI);
    PROTECT(NI);
    PROTECT(tolnI);
    PROTECT(debugI);
    PROTECT(rescale_variance_boolI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    int N = *(INTEGER(NI));
    int toln = *(INTEGER(tolnI));
    int debug = *INTEGER(debugI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);



    // SEXP Xsexp = PROTECT(allocVector(REALSXP, n*p));
    // double * X = REAL(Xsexp);

    // max vals for no partial sum
    SEXP maxvalSEXP = PROTECT(allocVector(REALSXP, N));
    double * maxval = REAL(maxvalSEXP);
    //memset(maxvals1, 0, sizeof(double)*len_as*N);

    for (int i = 0; i < N; ++i)
    {
        maxval[i] = NEGINF;
    }

    //memset(maxvals2, 0, sizeof(double)*len_as*N);

    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vectorSEXP = PROTECT(allocVector(REALSXP, maxlen));

   
    double * vector = REAL(vectorSEXP);
    memset(vector, 0, sizeof(double)*maxlen);

    SEXP XSEXP = PROTECT(allocVector(REALSXP, p * (n)));
    double *X = REAL(XSEXP); // p \times (n+1). first col is 0
    memset(X, 0, p*(n)*sizeof(double));
    

    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0

    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p);

    int s = -1;
    int e = n-1;
    double tmp = 0;
    double tmp2 = 0;
    int i= 0;
    int j= 0;
    for (int k = 0; k < N; ++k)
    {
        
        // generate X
        GetRNGstate();

        for (j = 0; j < n; ++j)
        {
            for (i = 0; i < p; ++i)
            {
                X[cord_spec(i,j,p)] = norm_rand();
            }
        }
        PutRNGstate();
        // normalize / rescale variance
        if (rescale_variance_bool)
        {
            rescale_variance(X,NULL, n, p, vector);
        }
        
        memset(cumsums, 0, p*sizeof(double));

        for (i = 0; i < p; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                cumsums[cord_spec(i,j+1,p)] = X[cord_spec(i,j,p)] +cumsums[cord_spec(i,j, p)];
            }
        }

        for (j = 0; j < e; ++j)
        {
            singleCUSUM(cumsums, cusum, s, e, p, j);
            //memset(tmpvec, 0, sizeof(double)*len_as);
            
            // first aggregate thresholded CUSUMs
            for (i = 0; i < p; ++i){
                tmp2 = cusum[cord_spec(i,0, p)]*cusum[cord_spec(i,0, p)];
                if (tmp2>maxval[k])
                {
                    maxval[k] = tmp2;
                }

            }
        }

            


    }


    //sort_k_largest(maxval, toln, 0, N);
    R_qsort(maxval, 1,N);
    SEXP maxSEXP = PROTECT(allocVector(REALSXP, 1));
    double * max = REAL(maxSEXP);
    //*max = maxval[toln-1];
    *max = maxval[N-toln];
    //Rprintf("%d\n", N-toln);



    UNPROTECT(12);
    return(maxSEXP);


}
