#include "header.h"

SEXP C_get_log_permanents(SEXP XSEXP, SEXP aSEXP, SEXP bSEXP, SEXP nSEXP, SEXP TSEXP, SEXP debugSEXP){

	PROTECT(XSEXP);
	PROTECT(aSEXP);
	PROTECT(bSEXP);
	PROTECT(nSEXP);
	PROTECT(TSEXP);
	PROTECT(debugSEXP);
	
	double * X = REAL(XSEXP);
	double * a = REAL(aSEXP);
	double * b = REAL(bSEXP);
	int n = *INTEGER(nSEXP);
	int T = *INTEGER(TSEXP);
	int debug = *INTEGER(debugSEXP);

	
	

	R_qsort(a, 1, n);
	R_qsort(b, 1, n);

	

	

	SEXP logpermsSEXP = PROTECT(allocVector(REALSXP, T));
	double * logperms = REAL(logpermsSEXP);
	memset(logperms, 0, sizeof(double)*T);


	SEXP a_union_bSEXP = PROTECT(allocVector(REALSXP, 2*n));
	double * a_union_b = REAL(a_union_bSEXP);
	memset(a_union_b, 0, sizeof(double)*2*n);

	int len_a_union_b =0;

	get_union(n, a, b, &len_a_union_b, a_union_b);

	SEXP alphaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP betaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP gammaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP logfactorialsSEXP = PROTECT(allocVector(REALSXP, n+1));
	SEXP mSEXP = PROTECT(allocVector(INTSXP, 1));
	SEXP k_SEXP = PROTECT(allocVector(INTSXP, 1));

	
	int * alpha = INTEGER(alphaSEXP);
	int * beta = INTEGER(betaSEXP);
	int * gamma = INTEGER(gammaSEXP);
	double * log_factorials = REAL(logfactorialsSEXP);
	int * m = INTEGER(mSEXP);
	int * k = INTEGER(k_SEXP);


	dictionary * new_log_subperms = init_dictionary(n);
	dictionary * old_log_subperms = init_dictionary(n);

	
	memset(alpha, 0, sizeof(int)*n);
	memset(beta, 0, sizeof(int)*n);
	memset(gamma, 0, sizeof(int)*n);
	memset(log_factorials, 0, sizeof(double)*(n+1));
	memset(m, 0, sizeof(int));
	memset(k, 0, sizeof(int));

	log_factorials[0]=0.0;
	for (int i = 1; i <= n; ++i)
	{
		log_factorials[i] = log_factorials[i-1] +log((double)(i));
	}

	
	
	SEXP historySEXP = PROTECT(allocVector(INTSXP, 3*n));
	SEXP amount_historySEXP = PROTECT(allocVector(INTSXP, 6*n));

	int * history = INTEGER(historySEXP);
	int * amount_history = INTEGER(amount_historySEXP);

	memset(history, 0, sizeof(int)*3*n);
	memset(amount_history, 0, sizeof(int)*6*n);

	for (int t = 0; t < T; ++t)
	{
		Rprintf("t = %d\n",t);
		double * x = X + (t*n);

		R_qsort(x, 1, n);

		if(!nonzero_perm(x, a,  b, n)){
			logperms[t] = -1;
			continue;
		}
		memset(alpha, 0, sizeof(int)*n);
		memset(beta, 0, sizeof(int)*n);
		memset(gamma, 0, sizeof(int)*n);
		memset(m, 0, sizeof(int));
		memset(k, 0, sizeof(int));

		get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
	    beta, gamma,  k, m, debug);


	    if(debug){
	    	Rprintf("len_a_union_b = %d\n", len_a_union_b);
	    	Rprintf("x:\n");
	    	print_float_vector(n,x);
	    	Rprintf("a:\n");
	    	print_float_vector(n,a);
	    	Rprintf("b:\n");
	    	print_float_vector(n,b);
	    	Rprintf("a_union_b:\n");
	    	print_float_vector(2*n,a_union_b);
	    	Rprintf("len a_union_b:%d\n", len_a_union_b);
	    	Rprintf("alpha:\n");
	    	print_int_vector(n,  alpha);
	    	Rprintf("beta:\n");
	    	print_int_vector(n,  beta);
	    	Rprintf("gamma:\n");
	    	print_int_vector(n,  gamma);
	    	Rprintf("m:%d\n", *m);
	    	Rprintf("k:%d\n", *k);
	    	
	    }

		int history_len = 0;

	
		memset(history, 0, sizeof(int)*3*n);
		memset(amount_history, 0, sizeof(int)*6*n);


		Rprintf("REDUCING NOW\n");
		//return(XSEXP);
		reduction(alpha,  beta,  gamma, m, n, k, history,
				   amount_history, &history_len, debug);

		Rprintf("history len = %d\n", history_len);

		Rprintf("REDUCED SUBPERMS\n");
		sparse_get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
						log_factorials, n,  m, k);
		//Rprintf("RESULT:\n");
		//print_matrix(n+1, n+1, new_log_subperms);

		dictionary * tmp  = old_log_subperms;
		old_log_subperms = new_log_subperms;
		new_log_subperms = tmp;




		Rprintf("==========\nReverse reduction:\n==========\n");
		//Rprintf("old = %d\n", old_log_subperms);
		//Rprintf("new = %d\n", new_log_subperms);
		dictionary * the_log_subperms = sparse_reverse_reduction(old_log_subperms, new_log_subperms, alpha,
						   beta,  gamma, m,  n, k,  history,
				           amount_history, &history_len, log_factorials);

		


		
		double logperm =  Csparse_log_sum_exp(the_log_subperms);
		logperms[t] = logperm;
		Rprintf("logperm = %f\n", logperm);



	}
	free_dictionary(new_log_subperms);
	free_dictionary(old_log_subperms);

	UNPROTECT(16);
	return(logpermsSEXP);

}


int nonzero_perm(double * x, double * a, double * b, int n){

	for (int i = 0; i < n; ++i)
	{
		if(x[i]< a[i] || x[i] > b[i]){
			return 0;
		}
	}
	return 1;



}