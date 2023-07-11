#include "header.h"



void get_union(int n, double * a, double * b, int * len_a_union_b, double * a_union_b){

    int counter =0;

    int a_counter = 0;
    int b_counter = 0;


    if(a[0]<b[0]){
        a_union_b[counter++] = a[0];
        a_counter++;
    }else{
        a_union_b[counter++] = b[0];
        b_counter++;
    }

    while(( a_counter <n )|| ( b_counter < n)){
        if(a_counter>=n){
            if(a_union_b[counter-1] < b[b_counter]){
                a_union_b[counter++] = b[b_counter];
            }
            b_counter++;

        }
        else if(b_counter>=n){
            if(a_union_b[counter-1] < a[b_counter]){
                a_union_b[counter++] = a[a_counter];
            }
            a_counter++;

        }
        else if (a[a_counter]<b[b_counter]){
            if(a_union_b[counter-1] < a[a_counter]){
                a_union_b[counter++] = a[a_counter];
            }
            a_counter++;
            
        }
        else if(b[b_counter]<a[a_counter]){
            if(a_union_b[counter-1] < b[b_counter])
            {
                a_union_b[counter++] = b[b_counter];
            }
            b_counter++;
        }

        else{
            if (a_union_b[counter-1] < b[b_counter])
            {
                a_union_b[counter++] = b[b_counter];
            }
            a_counter++;
            
        }
    }
    (*len_a_union_b) = counter;

}
void get_alphabetagamma(double * x, int n, double * a, double * b, double * a_union_b, int len_a_union_b, int * alpha, 
    int * beta, int * gamma, int * k, int * m, int debug){


    // initialize alpha beta gamma to zero
    memset(alpha, 0, sizeof(int)*n);
    memset(beta, 0, sizeof(int)*n);
    memset(gamma, 0, sizeof(int)*n);
    
    *m = n;
    memset(k, 0, sizeof(int));

    

    int alpha_counter = 0;
    int curr_x_ind = 0;

    int cumsum_beta = 0;
    int cumsum_beta_prev = 0;

    int cumsum_gamma = 0;
    int cumsum_gamma_next = 0;

    for (int i = 0; i < len_a_union_b-1; ++i)
    {
        if(debug){
            Rprintf("a_union_b[%d] = %f\n", i, a_union_b[i]);

        }
        for (int j = 0; j < n; ++j)
        {
            if ( (a_union_b[i] <= x[j]) && (x[j]<= a_union_b[i+1]))
            {
                alpha[alpha_counter]+=1;
                curr_x_ind = j;
            }
        }
        if(debug){
            Rprintf("alpha = %d\n", alpha[alpha_counter]);

        }
        


        if(alpha[alpha_counter]>0){

            if(alpha_counter>=1){
                cumsum_beta = 0;
                for (int z = 0; z < n; ++z)
                {
                    if(b[z] < x[curr_x_ind]){
                        cumsum_beta++;
                    }else{break;}
                }
                beta[alpha_counter-1] = cumsum_beta - cumsum_beta_prev;
            }

            cumsum_gamma_next = 0;
            for (int z = n-1; z >=0; --z)
                {
                    if(a[z] > x[curr_x_ind]){
                        cumsum_gamma_next++;
                    }else{break;}
                }

            if(alpha_counter>=1){
                gamma[alpha_counter-1] = cumsum_gamma - cumsum_gamma_next;
            }

            cumsum_gamma = cumsum_gamma_next;


            
            cumsum_beta_prev = cumsum_beta;
            alpha_counter++;

        }
        if(debug){
            Rprintf("alphacounter = %d\n", alpha_counter);

        }
    }

    *k = alpha_counter;

    return;
}


 SEXP C_test_get_alphabetagamma(){


    int n = 4;

    SEXP xSEXP = PROTECT(allocVector(REALSXP, n));
    SEXP aSEXP = PROTECT(allocVector(REALSXP, n));
    SEXP bSEXP = PROTECT(allocVector(REALSXP, n));
    SEXP aunionbSEXP = PROTECT(allocVector(REALSXP, 2*n));
    SEXP alphaSEXP = PROTECT(allocVector(INTSXP, n));
    SEXP betaSEXP = PROTECT(allocVector(INTSXP, n));
    SEXP gammaSEXP = PROTECT(allocVector(INTSXP, n));
    SEXP mSEXP = PROTECT(allocVector(INTSXP, 1));
    SEXP k_SEXP = PROTECT(allocVector(INTSXP, 1));

    
    int * alpha = INTEGER(alphaSEXP);
    int * beta = INTEGER(betaSEXP);
    int * gamma = INTEGER(gammaSEXP);
    int * m = INTEGER(mSEXP);
    int * k = INTEGER(k_SEXP);
    double * x = REAL(xSEXP);
    double * a = REAL(aSEXP);
    double * b = REAL(bSEXP);
    double * aunionb = REAL(aunionbSEXP);

    int debug = 1;

    
    memset(alpha, 0, sizeof(int)*n);
    memset(beta, 0, sizeof(int)*n);
    memset(gamma, 0, sizeof(int)*n);
    memset(m, 0, sizeof(int));
    memset(k, 0, sizeof(int));
    memset(x, 0, sizeof(double)*n);
    memset(a, 0, sizeof(double)*n);
    memset(b, 0, sizeof(double)*n);
    memset(aunionb, 0, sizeof(double)*2*n);

    
    a[0] = 0.0;
    a[1] = 1.0;
    a[2] = 3.0;
    a[3] = 5.0;

    b[0] = 2.0;
    b[1] = 4.0;
    b[2] = 6.0;
    b[3] = 7.0;

    aunionb[0] = 0.0;
    aunionb[1] = 1.0;
    aunionb[2] = 2.0;
    aunionb[3] = 3.0;
    aunionb[4] = 4.0;
    aunionb[5] = 5.0;
    aunionb[6] = 6.0;
    aunionb[7] = 7.0;
    int len_a_union_b = 8;
    x[0] = 1.5;
    x[1] = 3.5;
    x[2] = 3.5;
    x[3] = 6.5;

    if(debug){
        Rprintf("x = ");
        print_float_vector(n,  x);
        Rprintf("a = ");
        print_float_vector(n,  a);
        Rprintf("b = ");
        print_float_vector(n,  b);
        Rprintf("a_union_b = ");
        print_float_vector(len_a_union_b,  aunionb);
        Rprintf("here is my *C code*\n");

    }

    get_alphabetagamma(x,n, a, b, aunionb, len_a_union_b, alpha, 
    beta, gamma, k, m,1);

    
    if(debug){
        Rprintf("k = %d\n", *k);

        Rprintf("alpha = ");
        print_int_vector(n,  alpha);

        Rprintf("beta = ");
        print_int_vector(n,  beta);

        Rprintf("gamma = ");
        print_int_vector(n,  gamma);
        
    }


    UNPROTECT(9);



    return(mSEXP);

}



 SEXP C_test_get_alphabetagamma_input(SEXP xSEXP,SEXP aSEXP, SEXP bSEXP, SEXP nSEXP, SEXP debugSEXP){


    //int n = 4;


    /*SEXP xSEXP = PROTECT(allocVector(REALSXP, n));
    SEXP aSEXP = PROTECT(allocVector(REALSXP, n));
    SEXP bSEXP = PROTECT(allocVector(REALSXP, n));*/

    PROTECT(xSEXP);
    PROTECT(aSEXP);
    PROTECT(bSEXP);
    PROTECT(nSEXP);
    PROTECT(debugSEXP);
    int n = *(INTEGER(nSEXP));
    int debug = *(INTEGER(debugSEXP));

    SEXP aunionbSEXP = PROTECT(allocVector(REALSXP, 2*n));
    SEXP alphaSEXP = PROTECT(allocVector(INTSXP, n));
    SEXP betaSEXP = PROTECT(allocVector(INTSXP, n));
    SEXP gammaSEXP = PROTECT(allocVector(INTSXP, n));
    SEXP mSEXP = PROTECT(allocVector(INTSXP, 1));
    SEXP k_SEXP = PROTECT(allocVector(INTSXP, 1));

    
    int * alpha = INTEGER(alphaSEXP);
    int * beta = INTEGER(betaSEXP);
    int * gamma = INTEGER(gammaSEXP);
    int * m = INTEGER(mSEXP);
    int * k = INTEGER(k_SEXP);
    double * x = REAL(xSEXP);
    double * a = REAL(aSEXP);
    double * b = REAL(bSEXP);
    double * aunionb = REAL(aunionbSEXP);

    

    
    memset(alpha, 0, sizeof(int)*n);
    memset(beta, 0, sizeof(int)*n);
    memset(gamma, 0, sizeof(int)*n);
    memset(m, 0, sizeof(int));
    memset(k, 0, sizeof(int));
    memset(aunionb, 0, sizeof(double)*2*n);

    int len_a_union_b = 0;
    
    
    get_union(n, a, b, &len_a_union_b,  aunionb);

    if(debug){
        Rprintf("x = ");
        print_float_vector(n,  x);
        Rprintf("a = ");
        print_float_vector(n,  a);
        Rprintf("b = ");
        print_float_vector(n,  b);
        Rprintf("a_union_b = ");
        print_float_vector(len_a_union_b,  aunionb);
        Rprintf("here is my *C code*\n");

    }
    
    get_alphabetagamma(x,n, a, b, aunionb, len_a_union_b, alpha, 
    beta, gamma, k, m,debug);

    
    if(debug){
        Rprintf("k = %d\n", *k);

        Rprintf("alpha = ");
        print_int_vector(n,  alpha);

        Rprintf("beta = ");
        print_int_vector(n,  beta);

        Rprintf("gamma = ");
        print_int_vector(n,  gamma);
        
    }


    SEXP ret = PROTECT(allocVector(VECSXP, 5)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, alphaSEXP);
    SET_VECTOR_ELT(ret, 1, betaSEXP);
    SET_VECTOR_ELT(ret, 2, gammaSEXP);
    SET_VECTOR_ELT(ret, 3, mSEXP);
    SET_VECTOR_ELT(ret, 4, k_SEXP);
    

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 5));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("alpha"));
    SET_STRING_ELT(names, 1, mkChar("beta"));
    SET_STRING_ELT(names, 2, mkChar("gamma"));
    SET_STRING_ELT(names, 3, mkChar("m"));
    SET_STRING_ELT(names, 4, mkChar("k"));
    

    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(13);



    return(ret);

}