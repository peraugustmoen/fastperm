#include "header.h"


//' This function applies soft tresholding
void internal_soft_thresh(double * x, int len, double lambda){
    double tmp;
    for (int i = 0; i < len; ++i)
    {
        tmp = fabs(x[i]) - lambda;
        if(tmp<0){
            tmp=0.0;
        }
        else{
            if(x[i]<0){
                tmp = -tmp;
            }
        }
        x[i] = tmp;
    }
}

SEXP soft_thresh(SEXP xI, SEXP lenI, SEXP lambdaI){
    PROTECT(xI);
    PROTECT(lenI);
    PROTECT(lambdaI);


    int len = *INTEGER(lenI);
    double lambda = *REAL(lambdaI);
    double * x = REAL(xI);

    UNPROTECT(2);


    //SEXP out = PROTECT(allocVector(REALSXP, 1));

    internal_soft_thresh(x, len, lambda);

    UNPROTECT(1);
    return(xI);
}



double * internal_sparse_svd(double * Z, int r1, int c1, double lambda, double eps, int maxiter,
                        double * mhat, double * mhatprod, double * v, double * v2,int debug){

    memcpy(mhat, Z, r1*c1*sizeof(double));
    internal_soft_thresh(mhat, r1*c1, lambda);

    //SEXP ret;
    //SEXP tmpret;
    double * projection;
    double * tmproj;
    double summ=0;
    double sumsq=0;
    if(r1<c1){

        internal_matmultrightT(mhat, mhatprod,r1, c1);
        projection = internal_power_method(mhatprod, r1, eps, maxiter, v, v2,debug);
        if(projection==NULL){
          return NULL;
        }
    }
    else{

        //mhatprodSEXP = PROTECT(allocVector(REALSXP, c1*c1));
        //double *mhatprod = REAL(mhatprodSEXP);
        internal_matmultleftT(mhat, mhatprod,r1, c1);

        tmproj = internal_power_method(mhatprod, c1, eps, maxiter,v,v2,debug);
        if(tmproj==NULL){
          return NULL;
        }
        //PROTECT(tmpret);
        //tmpretarr = REAL(tmpret);
        //ret = PROTECT(allocVector(REALSXP, c1));
        //tmparr = REAL(ret);
        projection = v;
        if(tmproj ==v){
            projection = v2;
        }
        internal_matmult(mhat, tmproj, projection, r1, c1, r1, 1);

        for (int i = 0; i < r1; ++i)
        {
            summ+=projection[i]*projection[i];
        }
        sumsq=sqrt(summ);
        for (int i = 0; i < r1; ++i)
        {
            projection[i] = projection[i]/sumsq;
        }


    }

    /*UNPROTECT(3);
    if (r1>=c1)
    {
        UNPROTECT(1);
    }*/
    return projection;
}

SEXP sparse_svd(SEXP ZI, SEXP r1I, SEXP c1I, SEXP lambdaI, SEXP epsI, SEXP maxiterI){
    PROTECT(ZI);
    PROTECT(r1I);
    PROTECT(c1I);
    PROTECT(lambdaI);
    PROTECT(epsI);
    PROTECT(maxiterI);

    double *Z= REAL(ZI);
    int r1 = *(INTEGER(r1I));
    int c1 = *(INTEGER(c1I));
    double lambda= *(REAL(lambdaI));
    double eps = *(REAL(epsI));
    int maxiter= *(INTEGER(maxiterI));
    UNPROTECT(5);

    int maxlen = r1;
    int minlen = c1;
    if(c1>r1){
        maxlen = c1;
        minlen = r1;
    }
    SEXP vec1SEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP vec2SEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP mhatSEXP = PROTECT(allocVector(REALSXP, r1*c1));
    SEXP mhatprodSEXP = PROTECT(allocVector(REALSXP, minlen*minlen));

    double * vec1 = REAL(vec1SEXP);
    double * vec2 = REAL(vec2SEXP);
    double * mhat = REAL(mhatSEXP);
    double * mhatprod = REAL(mhatprodSEXP);

    //internal_sparse_svd(double * Z, int r1, int c1, double lambda, double eps, int maxiter,
                        //double * mhat, double * mhatprod, double * v, double * v2)

    double* retarr = internal_sparse_svd(Z, r1, c1, lambda, eps, maxiter,
                        mhat, mhatprod, vec1, vec2,0);
    SEXP ret = vec2SEXP;
    if(retarr==vec1){
        ret = vec1SEXP;
    }
    UNPROTECT(5);
    return ret;

}





void internal_inspectOnSegment(double * cumsums, double * cusum, int * maxpos, double * maximum, int s, int e, int p, double lambda,
    double eps, int maxiter, double * mhat, double * mhatprod, double* v, double* v2,int debug){
    *maxpos  = e;
    *maximum = 0.0;
    if(e-s<2){
        return;
    }

    // compute CUSUM
    CUSUM(cumsums, cusum, s, e, p);

    // find sparse SVD
    //internal_sparse_svd(double * Z, int r1, int c1, double lambda, double eps, int maxiter,
                        //double * mhat, double * mhatprod, double * v, double * v2)
    double * projvec = internal_sparse_svd(cusum, p, e-s-1, lambda, eps, maxiter,
                        mhat, mhatprod, v, v2,debug);
    if(projvec==NULL){
      if(debug){
          Rprintf("inspecting segment, s=%d, e=%d resulted in NULL projection. lambda = %f.\n", s,e,lambda);
      }

      return;
    }
    double * projected = v;
    if(projvec==v){
        projected = v2;
    }
    //double * v = REAL(retval);
    int n = e-s;
    double tmp;
    internal_matmult(projvec,cusum,  projected,1, p, p, n-1);

    for (int i = 0; i < n-1; ++i)
    {
        //t = s+i+1;
        tmp = fabs(projected[i]);
        if(tmp> *maximum){
            *maximum = tmp;
            *maxpos = s+i+1;
        }
    }
    if(debug){
        Rprintf("inspecting segment, s=%d, e=%d, max_cusum = %f\n", s,e,*maximum);
    }


    return;
}

void cInspect_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter,
                double * maxval, double xi, double *cumsums, int* lens, int lenLens, double lambda,
                double eps, int maxiter, int * segstarts, double * maxcusums, int* maxpos, int K, double * cusum, double * mhat,
                double * mhatprod, double * v, double * v2, int debug,int * coordchg){
    if(debug){
        Rprintf("cInspectCall! s=%d, e=%d\n", s, e);
    }

    if(e-s<2*lens[0]){
        //Rprintf("segment too short\n");
        return;
    }
    int argmax = s;
    double maximum = 0;

    //int tmpargmax;
    //double tmpmaximum;
    double tmp;
    int len;
    int jump;
    //int found=0;

    int i;

    int j_max=0;
    int k_max = 0;
    for (int j = 0; j < lenLens; ++j)
    {
        len = lens[j];
        if(debug){
            Rprintf("j=%d, len = %d\n", j, len);
        }

        jump = len /K;
        if(jump==0){
            jump = 1;
        }

        if(e-s<2*len){
            break;
        }

        for (int k = 0; k < n; ++k)
        {
            i = segstarts[cord_spec(k,j,n)];

            if(i>e-2*len || i<-1){
                break;
            }
            else if (i<s)
            {
                continue;
            }

            if(debug){
                Rprintf("maxcusums[%d, %d] = %f\n", k , j , maxcusums[cord_spec(k,j,n)]);
            }

            if(maxcusums[cord_spec(k,j,n)]<=0.0){
                //this segment not computed!
                //inspectOnSegment(cumsums, cusum, &tmpargmax, &tmpmaximum, s, e, p, lambda,
                 //   eps, maxiter, projvec, cusum_proj);
                 //double * cumsums, double * cusum, int * maxpos, double * maximum,
                 internal_inspectOnSegment(cumsums, cusum, &(maxpos[cord_spec(k,j,n)]), &(maxcusums[cord_spec(k,j,n)]), i, i+2*len, p,
                 lambda,
                        eps, maxiter, mhat, mhatprod, v, v2,debug);
            }

            tmp = maxcusums[cord_spec(k,j,n)];
            if(tmp>maximum){
                maximum = tmp;
                argmax = maxpos[cord_spec(k,j,n)];
                j_max = j;
                k_max=k;
                //found=1;
            }


        }

        if(maximum>xi){
            break;
        }
    }
    if(debug){
        Rprintf("maximum=%f\n", maximum);
    }

    if(maximum > xi){
        if(debug){
            Rprintf("!!!!!! declared change-point in %d. val = %f, thresh =%f\n", argmax, maximum, xi);
        }
        // identify in which coordinates the change happens:
        i = segstarts[cord_spec(k_max,j_max,n)];
        len = lens[j_max];
        int ss = i;
        int ee = i+2*len;
        CUSUM(cumsums, cusum, ss, ee, p);
        double * projvec = internal_sparse_svd(cusum, p, ee-ss-1, lambda, eps, maxiter,
                        mhat, mhatprod, v, v2,debug);
        for (int zz = 0; zz < p; ++zz)
        {
            if(fabs(projvec[zz])>1e-6){
                coordchg[cord_spec(zz,*changepoint_counter_ptr, p)]=1;
            }
        }
        
        changepoints[*changepoint_counter_ptr] = argmax;
        depthcounter[*changepoint_counter_ptr] = depth;
        maxval[*changepoint_counter_ptr] = maximum;
        (*changepoint_counter_ptr)++;
        //cInspect_call(x, s, argmax, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,  maxval, threshold,
        //        adaptTresh, cumsums, lens, lenLens, K, cusum,  mhat, mhatprod, v, v2);
        cInspect_call(x, s, argmax, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval,xi, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxcusums, maxpos, K,  cusum, mhat,
                mhatprod, v, v2,debug,coordchg);
        //cInspect_call(x, argmax, e, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,  maxval, threshold,
        //        adaptTresh, cumsums, lens, lenLens, K, cusum, mhat, mhatprod, v, v2);
        cInspect_call(x, argmax, e, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval,xi, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxcusums, maxpos, K,  cusum, mhat,
                mhatprod, v, v2,debug,coordchg);

    }

    return;
}

SEXP cInspect(SEXP XI,SEXP nI, SEXP pI,SEXP xiI, SEXP lensI,SEXP lenLensI,SEXP KI,
    SEXP epsI, SEXP lambdaI, SEXP maxiterI, SEXP rescale_variance_boolI, SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(lensI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(xiI);
    PROTECT(lenLensI);
    PROTECT(KI);
    PROTECT(epsI);
    PROTECT(lambdaI);
    PROTECT(maxiterI);
    PROTECT(rescale_variance_boolI);
    PROTECT(debugI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double xi = *(REAL(xiI));
    int *lens = INTEGER(lensI); /////// dobbeltsjekk at int er like stor som INTEGER!!!!!!!
    int lenLens = *(INTEGER(lenLensI));
    int K = *(INTEGER(KI));
    double eps = *REAL(epsI);
    double lambda = *REAL(lambdaI);
    int maxiter = *INTEGER(maxiterI);
    int debug = *INTEGER(debugI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);
    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("lambda = %f\n", lambda);
    }

    SEXP out = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out); //pointer to array
    memset(changepoints, -1, sizeof(int)*n);
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    SEXP maxvalSEXP = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(maxvalSEXP); //pointer to array
    memset(maxval, 0, sizeof(double)*n);
    SEXP depthcounterSEXP= PROTECT(allocVector(INTSXP, n));
    int * depthcounter = INTEGER(depthcounterSEXP); //pointer to array
    memset(depthcounter, 0, sizeof(int)*n);
    SEXP coordschgSEXP = PROTECT(allocVector(INTSXP, n*p));
    int * coordchg = INTEGER(coordschgSEXP);
    memset(coordchg, 0, sizeof(int)*n*p);
    // first we compute all the cumulative sums of all
    // coordinates:

    

    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p*(n)));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p*(n));
    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vSEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP v2SEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP mhatSEXP = PROTECT(allocVector(REALSXP, p*n));
    SEXP mhatprodSEXP = PROTECT(allocVector(REALSXP, minlen*minlen));

    double * v = REAL(vSEXP);
    memset(v, 0, sizeof(double)*maxlen);
    double * v2 = REAL(v2SEXP);
    memset(v2, 0, sizeof(double)*maxlen);
    double * mhat = REAL(mhatSEXP);
    memset(mhat, 0, sizeof(double)*p*n);
    double * mhatprod = REAL(mhatprodSEXP);
    memset(mhatprod, 0, sizeof(double)*minlen*minlen);

    SEXP scalesSEXP = PROTECT(allocVector(REALSXP, p));
    double * scales = REAL(scalesSEXP);
    for (int i = 0; i < p; ++i)
    {
        scales[i] = 1;
    }


    if (rescale_variance_bool)
    {
        rescale_variance(X, scales, n, p, v);
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


    // record all segment starts
    // we don't compute all cusums here but do it on the fly so we dont have to do anything
    // unecessary


    SEXP maxcusumsSEXP = PROTECT(allocVector(REALSXP, n*lenLens));
    double * maxcusums = REAL(maxcusumsSEXP);
    memset(maxcusums, 0, sizeof(double)*n*lenLens);
    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxpos = INTEGER(maxposSEXP);
    memset(maxpos, 0, sizeof(int)*n*lenLens);
    SEXP segstartsSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * segstarts = INTEGER(segstartsSEXP);
    memset(segstarts, -2, sizeof(int)*n*lenLens);

    int len;
    int jump;
    int counter = 0;
    for (int j = 0; j < lenLens; ++j)
    {
        counter=0;

        len = lens[j]; //len here is -1 of the len in text
        jump = len/K;
        if(jump<1){
            jump=1;
        }

        for (int i = -1; i < (n-2*len); i+=jump)
        {
            //cord_spec(r,c, D) ((r) + (D)*(c))
            //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxcusums[cord_spec(i,j,n)]));
            segstarts[cord_spec(counter++,j,n)] = i;
            if(debug){
                Rprintf("segstarts[%d, %d] = %d\n",counter-1, j, i );
            }


        }
    }



    cInspect_call(X, -1, n-1, n, p, 1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval, xi, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxcusums, maxpos, K, cusum, mhat,
                mhatprod, v, v2,debug,coordchg);

    // return:
    SEXP changepointnumSEXP = PROTECT(allocVector(INTSXP, 1));
    int * changepointnum = INTEGER(changepointnumSEXP);
    *changepointnum = changepoint_counter;

    SEXP ret = PROTECT(allocVector(VECSXP, 6)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out);
    SET_VECTOR_ELT(ret, 1, maxvalSEXP);
    SET_VECTOR_ELT(ret, 2, depthcounterSEXP);
    SET_VECTOR_ELT(ret, 3, coordschgSEXP);
    SET_VECTOR_ELT(ret, 4, changepointnumSEXP);
    SET_VECTOR_ELT(ret, 5, scalesSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 6));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("CUSUMval"));
    SET_STRING_ELT(names, 2, mkChar("depth"));
    SET_STRING_ELT(names, 3, mkChar("coordinate"));
    SET_STRING_ELT(names, 4, mkChar("changepointnumber"));
    SET_STRING_ELT(names, 5, mkChar("scales"));

    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(29);
    return ret;
}


SEXP cInspect_single(SEXP XI,SEXP nI, SEXP pI,SEXP xiI, 
    SEXP epsI, SEXP lambdaI, SEXP maxiterI, SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(xiI);
   
    PROTECT(epsI);
    PROTECT(lambdaI);
    PROTECT(maxiterI);
    PROTECT(debugI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double xi = *(REAL(xiI));
    double eps = *REAL(epsI);
    double lambda = *REAL(lambdaI);
    int maxiter = *INTEGER(maxiterI);
    int debug = *INTEGER(debugI);
    UNPROTECT(7); // unprotecting all except X and lens
    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("lambda = %f\n", lambda);
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
    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vSEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP v2SEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP mhatSEXP = PROTECT(allocVector(REALSXP, p*n));
    SEXP mhatprodSEXP = PROTECT(allocVector(REALSXP, minlen*minlen));

    double * v = REAL(vSEXP);
    memset(v, 0, sizeof(double)*maxlen);
    double * v2 = REAL(v2SEXP);
    memset(v2, 0, sizeof(double)*maxlen);
    double * mhat = REAL(mhatSEXP);
    memset(mhat, 0, sizeof(double)*p*n);
    double * mhatprod = REAL(mhatprodSEXP);
    memset(mhatprod, 0, sizeof(double)*minlen*minlen);



    // record all segment starts
    // we don't compute all cusums here but do it on the fly so we dont have to do anything
    // unecessary

    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, 1));
    int * maxpos = INTEGER(maxposSEXP);
    *maxpos = 0;

    SEXP maximumSEXP = PROTECT(allocVector(REALSXP, 1));
    double * maximum = REAL(maximumSEXP);
    *maximum = -1000000000000000000000.0;
    int s = -1;
    int e = n-1;

    internal_inspectOnSegment(cumsums, cusum, maxpos, maximum, s, e, p,
                 lambda,eps, maxiter, mhat, mhatprod, v, v2,debug);

    SEXP ret = PROTECT(allocVector(VECSXP, 2)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, maxposSEXP);
    SET_VECTOR_ELT(ret, 1, maximumSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("pos"));
    SET_STRING_ELT(names, 1, mkChar("cusumval"));


    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(11);
    return ret;
}


double internal_Inspect_single(double * X,int n, int p, 
    double eps, double lambda, int maxiter, int debug){
    // X : p \times n
    




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
    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vSEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP v2SEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP mhatSEXP = PROTECT(allocVector(REALSXP, p*n));
    SEXP mhatprodSEXP = PROTECT(allocVector(REALSXP, minlen*minlen));

    double * v = REAL(vSEXP);
    memset(v, 0, sizeof(double)*maxlen);
    double * v2 = REAL(v2SEXP);
    memset(v2, 0, sizeof(double)*maxlen);
    double * mhat = REAL(mhatSEXP);
    memset(mhat, 0, sizeof(double)*p*n);
    double * mhatprod = REAL(mhatprodSEXP);
    memset(mhatprod, 0, sizeof(double)*minlen*minlen);



    // record all segment starts
    // we don't compute all cusums here but do it on the fly so we dont have to do anything
    // unecessary

    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, 1));
    int * maxpos = INTEGER(maxposSEXP);
    *maxpos = 0;


    double maximum = -1000000000000000000000.0;
    int s = -1;
    int e = n-1;

    internal_inspectOnSegment(cumsums, cusum, maxpos, &maximum, s, e, p,
                 lambda,eps, maxiter, mhat, mhatprod, v, v2,debug);

    UNPROTECT(7);
    return maximum;
}

/*SEXP cInspect_test_calibrate(SEXP nI, SEXP pI, SEXP NI, SEXP tolnI, SEXP lambdaI, SEXP epsI,
    SEXP maxiterI, SEXP rescale_variance_boolI, SEXP debugI){

    PROTECT(nI);
    PROTECT(pI);
    PROTECT(NI);
    PROTECT(tolnI);
    PROTECT(lambdaI);
    PROTECT(epsI);
    PROTECT(maxiterI);
    PROTECT(debugI);
    PROTECT(rescale_variance_boolI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    int N = *(INTEGER(NI));
    int toln = *(INTEGER(tolnI));
    double lambda = *(REAL(lambdaI));
    double eps = *(REAL(epsI));
    int maxiter = *(INTEGER(maxiterI));
    int debug = *INTEGER(debugI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);



    if(debug){
        Rprintf("n = %d\np = %d\nlambda = %f\ntoln = %d\nrescale = %d\n",n,p,lambda,toln,rescale_variance_bool);
    }
    // SEXP Xsexp = PROTECT(allocVector(REALSXP, n*p));
    // double * X = REAL(Xsexp);

    // max vals for no partial sum
    SEXP maxvals1SEXP = PROTECT(allocVector(REALSXP, N));
    double * maxvals1 = REAL(maxvals1SEXP);
    memset(maxvals1, 0, sizeof(double)*N);

   

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
    double *X = REAL(XSEXP);
    memset(X, 0, sizeof(double)*p*n);
  
    double maxx = -1000000000.0;
    GetRNGstate();
    int i=0;
    int j = 0;

    for (int k = 0; k < N; ++k)
    {
        
        // generate X

        for (j = 0; j < n; ++j)
        {
            for (i = 0; i < p; ++i)
            {
                X[cord_spec(i,j,p)] = norm_rand();
            }
        }

        // normalize / rescale variance
        if (rescale_variance_bool)
        {
            rescale_variance(X,NULL, n, p, vector);
        }
        
        maxvals1[k] = internal_Inspect_single(X,n, p, 
                        eps, lambda, maxiter, debug);

        if(maxvals1[k]>maxx){
            maxx = maxvals1[k];
        }

        if(debug){
            Rprintf("maxval[%d] = %f. cur max = %f\n", k , maxvals1[k],maxx);
        }



    }

    PutRNGstate();
    SEXP thresholds1SEXP = PROTECT(allocVector(REALSXP, 1));
    double * thresholds1 = REAL(thresholds1SEXP);


    sort_k_largest(maxvals1 , toln, 0, N);
        
    thresholds1[0] = maxvals1[cord_spec(toln-1, 0, N)];


    SEXP ret = PROTECT(allocVector(VECSXP, 1)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, thresholds1SEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 1));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("xi"));

    setAttrib(ret, R_NamesSymbol, names);


    UNPROTECT(15);
    return(ret);


}*/



SEXP cInspect_calibrate(SEXP nI, SEXP pI, SEXP NI, SEXP tolnI,SEXP lensI,SEXP lenLensI,SEXP KI, SEXP lambdaI,
                        SEXP epsI,SEXP maxiterI, SEXP rescale_variance_boolI, SEXP debugI){
    // X : p \times n
    PROTECT(lensI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(NI);
    PROTECT(tolnI);
    PROTECT(lenLensI);
    PROTECT(KI);
    PROTECT(debugI);
    PROTECT(lambdaI);
    PROTECT(maxiterI);
    PROTECT(rescale_variance_boolI);
    PROTECT(epsI);

    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    int N = *(INTEGER(NI));
    int toln = *(INTEGER(tolnI));
    int *lens = INTEGER(lensI); /////// dobbeltsjekk at int er like stor som INTEGER!!!!!!!
    int lenLens = *(INTEGER(lenLensI));
    int K = *(INTEGER(KI));
    double lambda = *(REAL(lambdaI));
    int debug = *INTEGER(debugI);
    int maxiter = *INTEGER(maxiterI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);
    double eps = *REAL(epsI);

    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("n = %d\n", n);
    }

    // max val for each iter
    SEXP maxvalsSEXP = PROTECT(allocVector(REALSXP, N));
    double * maxvals = REAL(maxvalsSEXP);
    //memset(maxvals1, 0, sizeof(double)*len_as*N);

    

    for (int i = 0; i < N; ++i)
    {
        maxvals[i] = NEGINF;
    }



    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));


    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p*(n)));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p*(n));
    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vSEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP v2SEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP mhatSEXP = PROTECT(allocVector(REALSXP, p*n));
    SEXP mhatprodSEXP = PROTECT(allocVector(REALSXP, minlen*minlen));

    double * v = REAL(vSEXP);
    memset(v, 0, sizeof(double)*maxlen);
    double * v2 = REAL(v2SEXP);
    memset(v2, 0, sizeof(double)*maxlen);
    double * mhat = REAL(mhatSEXP);
    memset(mhat, 0, sizeof(double)*p*n);
    double * mhatprod = REAL(mhatprodSEXP);
    memset(mhatprod, 0, sizeof(double)*minlen*minlen);

    SEXP segstartsSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * segstarts = INTEGER(segstartsSEXP);
    //memset(segstarts, -2, sizeof(int)*n*lenLens);
    for (size_t i = 0; i < n*lenLens; i++) {
      segstarts[i] = -2;
    }
    int len;
    int jump;
    int counter = 0;

    for (int j = 0; j < lenLens; ++j)
    {
        counter=0;

        len = lens[j]; //len here is -1 of the len in text
        jump = len/K;
        if(jump<1){
            jump=1;
        }

        for (int i = -1; i < (n-2*len); i+=jump)
        {
            //cord_spec(r,c, D) ((r) + (D)*(c))
            //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxvalues[cord_spec(i,j,n)]));
            segstarts[cord_spec(counter++,j,n)] = i;
            if(debug){
                //Rprintf("segstarts[%d, %d] = %d\n",counter-1, j, i );
            }


        }
    }

    int s = -1;
    int e = n-1;
    int ss;
    int ee;
    int ll;
    int hh;
    int i;
    int j;
    int k;

    double tmpmax = 0.0;
    int pos = 0;

    /*if(debug){
        return NULL;
    }*/
    GetRNGstate();

    SEXP XSEXP = PROTECT(allocVector(REALSXP, p * (n)));
    double *X = REAL(XSEXP); // p \times (n+1). first col is 0
    memset(X, 0, p*(n)*sizeof(double));

    for (k = 0; k < N; ++k)
    {

        // generate X

        for (j = 0; j < n; ++j)
        {
            for (i = 0; i < p; ++i)
            {
                X[cord_spec(i,j,p)] = norm_rand();
            }
        }

        // normalize / rescale variance
        if (rescale_variance_bool)
        {
            rescale_variance(X, NULL, n, p, v);
        }
        
        memset(cumsums, 0, p*sizeof(double));

        for (i = 0; i < p; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                cumsums[cord_spec(i,j+1,p)] = X[cord_spec(i,j,p)] +cumsums[cord_spec(i,j, p)];
            }
        }


        for ( ll = 0; ll < lenLens; ++ll)
        {
            len = lens[ll];
            if(debug){
                Rprintf("ll=%d, len = %d\n", ll, len);
            }



            if(e-s<2*len){
                break;
            }

            for (hh = 0; hh < n; ++hh)
            {
                ss = segstarts[cord_spec(hh,ll,n)];
                ee = ss+2*len;
                if(debug){
                    Rprintf("i= %d\n", hh);
                }
                if( ee >n || ss<-1){
                    if(debug){
                        Rprintf("i= %d is shhipped\n", hh);
                    }
                    break;
                }





                internal_inspectOnSegment(cumsums, cusum, &pos, &tmpmax, ss, ee, p, lambda,
                    eps, maxiter, mhat, mhatprod, v, v2, debug);

                if(tmpmax >maxvals[k]){
                    maxvals[k] = tmpmax;
                }
                
            }
        }
    }



    PutRNGstate();
    SEXP thresholdSEXP = PROTECT(allocVector(REALSXP, 1 ));
    double * threshold = REAL(thresholdSEXP);

        sort_k_largest(maxvals , toln, 0, N);
        *threshold = maxvals[cord_spec(toln-1, 0, N)];

  

    SEXP ret = PROTECT(allocVector(VECSXP, 1)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, thresholdSEXP);


    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 1));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("max_value"));

    setAttrib(ret, R_NamesSymbol, names);


    UNPROTECT(24);
    return(ret);



}





SEXP cInspect_test_calibrate(SEXP nI, SEXP pI, SEXP NI, SEXP tolnI, SEXP lambdaI,
                        SEXP epsI,SEXP maxiterI, SEXP rescale_variance_boolI, SEXP debugI){
    // X : p \times n
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(NI);
    PROTECT(tolnI);
    PROTECT(debugI);
    PROTECT(lambdaI);
    PROTECT(maxiterI);
    PROTECT(rescale_variance_boolI);
    PROTECT(epsI);

    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    int N = *(INTEGER(NI));
    int toln = *(INTEGER(tolnI));
    double lambda = *(REAL(lambdaI));
    int debug = *INTEGER(debugI);
    int maxiter = *INTEGER(maxiterI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);
    double eps = *REAL(epsI);

    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("n = %d\n", n);
    }

    /*Rprintf("toln = %d\n",toln);
    Rprintf("rescale = %d\n",rescale_variance_bool);*/

    // max val for each iter
    SEXP maxvalsSEXP = PROTECT(allocVector(REALSXP, N));
    double * maxvals = REAL(maxvalsSEXP);
    //memset(maxvals1, 0, sizeof(double)*len_as*N);

    

    for (int i = 0; i < N; ++i)
    {
        maxvals[i] = NEGINF;
    }



    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));


    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p*(n)));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p*(n));
    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vSEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP v2SEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP mhatSEXP = PROTECT(allocVector(REALSXP, p*n));
    SEXP mhatprodSEXP = PROTECT(allocVector(REALSXP, minlen*minlen));

    double * v = REAL(vSEXP);
    memset(v, 0, sizeof(double)*maxlen);
    double * v2 = REAL(v2SEXP);
    memset(v2, 0, sizeof(double)*maxlen);
    double * mhat = REAL(mhatSEXP);
    memset(mhat, 0, sizeof(double)*p*n);
    double * mhatprod = REAL(mhatprodSEXP);
    memset(mhatprod, 0, sizeof(double)*minlen*minlen);

    /*SEXP segstartsSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * segstarts = INTEGER(segstartsSEXP);
    //memset(segstarts, -2, sizeof(int)*n*lenLens);
    for (size_t i = 0; i < n*lenLens; i++) {
      segstarts[i] = -2;
    }
    int len;
    int jump;
    int counter = 0;

    for (int j = 0; j < lenLens; ++j)
    {
        counter=0;

        len = lens[j]; //len here is -1 of the len in text
        jump = len/K;
        if(jump<1){
            jump=1;
        }

        for (int i = -1; i < (n-2*len); i+=jump)
        {
            //cord_spec(r,c, D) ((r) + (D)*(c))
            //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxvalues[cord_spec(i,j,n)]));
            segstarts[cord_spec(counter++,j,n)] = i;
            if(debug){
                //Rprintf("segstarts[%d, %d] = %d\n",counter-1, j, i );
            }


        }
    }*/

    int s = -1;
    int e = n-1;
    int ss;
    int ee;
    int ll;
    int hh;
    int i;
    int j;
    int k;

    double tmpmax = 0.0;
    int pos = 0;

    /*if(debug){
        return NULL;
    }*/
    

    SEXP XSEXP = PROTECT(allocVector(REALSXP, p * (n)));
    double *X = REAL(XSEXP); // p \times (n+1). first col is 0
    memset(X, 0, p*(n)*sizeof(double));

    for (k = 0; k < N; ++k)
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
         

        //Rprintf("XXX = %f\n", X[n*p-1]);
        // normalize / rescale variance
        if (rescale_variance_bool)
        {
            rescale_variance(X, NULL, n, p, v);
        }

        //Rprintf("XXX rescaled = %f\n", X[n*p-1]);
        
        memset(cumsums, 0, p*sizeof(double));

        for (i = 0; i < p; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                cumsums[cord_spec(i,j+1,p)] = X[cord_spec(i,j,p)] +cumsums[cord_spec(i,j, p)];
            }
        }



        ss = -1;
        ee = n-1;
        /*if(debug){
            Rprintf("i= %d\n", hh);
        }*/






        internal_inspectOnSegment(cumsums, cusum, &pos, &tmpmax, ss, ee, p, lambda,
            eps, maxiter, mhat, mhatprod, v, v2, debug);



        if(tmpmax >maxvals[k]){
            maxvals[k] = tmpmax;
        }
          
            
    }



    
    SEXP thresholdSEXP = PROTECT(allocVector(REALSXP, 1 ));
    double * threshold = REAL(thresholdSEXP);

    R_qsort(maxvals , 1,N);

    //Rprintf("%d\n",toln);
    //Rprintf("%d\n",N-toln);
    *threshold = maxvals[cord_spec(N-toln, 0, N)];
/*    sort_k_largest(maxvals , toln, 0, N);
    *threshold = maxvals[cord_spec(toln-1, 0, N)];
*/
  

    SEXP ret = PROTECT(allocVector(VECSXP, 1)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, thresholdSEXP);


    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 1));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("max_value"));

    setAttrib(ret, R_NamesSymbol, names);


    UNPROTECT(20);
    return(ret);



}



