#include "header.h"




// reverse top trim
void sparse_reverse_tt(dictionary * old_log_subperms,  dictionary * new_log_subperms,
						int amount, int amount2,
						int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m,
						int * k){

	int old_r=0;
	int old_s=0;
	int new_r = 0;
	int new_s = 0;

	double log_lower = 0;
	double log_upper = 0;
	double new_value = 0;
	//double log_lower = 0;

	nullset_dictionary(new_log_subperms);

	for (int z = 0; z < (*old_log_subperms).used_len; ++z)
	{
		pair old_pair = (*old_log_subperms).array[z];

		old_r = old_pair.x;
		old_s = old_pair.y;

		new_r = old_r + amount;
		new_s = old_s;



		// check bounds
		if(new_r< amount || new_r> alpha[0]){
			continue;
		}
		if(new_s<0 || new_s > alpha[(*k)-1]){
			continue;
		}

		pair new_pair =  { .x=new_r, .y = new_s};




		log_lower = log_factorials[alpha[0] - new_r];
		log_upper = log_factorials[amount + alpha[0] - new_r];

		new_value = log_upper - log_lower + (*old_log_subperms).value_array[z];

		
		add_to_dictionary(new_log_subperms, new_pair, new_value);
	

	}

	// change alpha beta osv osv

	beta[0] = amount;

	*m = (*m) + amount;

	return;

}



// Reverse bottom trim
void sparse_reverse_bt(dictionary * old_log_subperms,  dictionary * new_log_subperms,
						int amount, int amount2,
						int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m,
						int * k){

	int old_r=0;
	int old_s=0;
	int new_r = 0;
	int new_s = 0;

	double log_lower = 0;
	double log_upper = 0;
	double new_value = 0;
	//double log_lower = 0;

	nullset_dictionary(new_log_subperms);

	for (int z = 0; z < (*old_log_subperms).used_len; ++z)
	{
		pair old_pair = (*old_log_subperms).array[z];

		old_r = old_pair.x;
		old_s = old_pair.y;

		new_r = old_r;
		new_s = old_s +amount;



		// check bounds
		if(new_r< 0|| new_r> alpha[0]){
			continue;
		}
		if(new_s< amount || new_s > alpha[(*k)-1]){
			continue;
		}

		pair new_pair =  { .x=new_r, .y = new_s};




		log_upper = log_factorials[ alpha[(*k)-1] +amount - new_s];

		log_lower = log_factorials[alpha[(*k)-1] - new_s];

		new_value = log_upper - log_lower + (*old_log_subperms).value_array[z];

		
		add_to_dictionary(new_log_subperms, new_pair, new_value);
	

	}

	// change alpha beta osv osv

	gamma[(*k) - 2] = amount;
	*m = (*m) + amount;

	return;

}

// Reverse left merge
void sparse_reverse_lm(dictionary * old_log_subperms,  dictionary * new_log_subperms,
						int amount, int amount2,
						int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m,
						int * k){


	int s = 0;
	int l=0;

	double log_upper1=0;
	double log_upper2=0;
	double log_upper3=0;
	double log_lower1=0;
	double log_lower2=0;
	double log_lower3=0;
	double log_upper4=0;
	double log_lower4a=0;
	double log_lower4b=0;

	double value = 0;
	double new_value = 0;
	//double log_lower = 0;
	clock_t start, end;
	pair new_pair = {.x = 0,.y=0};

	nullset_dictionary(new_log_subperms);

	for (int z = 0; z < (*old_log_subperms).used_len; ++z)
	{
		
		
		pair old_pair = (*old_log_subperms).array[z];

		l = old_pair.x;
		s = old_pair.y;

		value = (*old_log_subperms).value_array[z];
		// check bounds
		// 0<=r<=amount
		// 0<= s <= alpha[(*k)-1]
		// r<=l<=r + amount2
		
		new_pair.y = s;

		for (int r = MAX(l-amount2, 0); r <= MIN(l, amount); ++r)
		{	
			/*if(z<10){
				start = clock();
			}*/
			log_upper1 = log_factorials[amount2];
			log_upper2 = log_factorials[amount];
			log_upper3 = log_factorials[amount + amount2];
			log_upper4 = log_factorials[l];
			log_lower1 = log_factorials[amount2 - l +r];
			log_lower2 = log_factorials[amount - r];
			log_lower3 = log_factorials[ amount +amount2 -l];
			log_lower4a = log_factorials[l-r];
			log_lower4b = log_factorials[r];
			/*if(z<10){
				end = clock();
				Rprintf("used %f clock cycles for computing loglowerupperosv\n", ((double) (end - start)));
			}*/	
			new_value = log_upper1 + log_upper2 - log_upper3 + log_upper4
 								- log_lower1 - log_lower2 + log_lower3
 								- log_lower4a -log_lower4b + value;
 			new_pair.x = r;
 			/*if(z<10){
				start = clock();
			}*/
 			update_dict(new_pair, new_value, new_log_subperms);
 			/*if(z<10){
				end = clock();
				Rprintf("used %f clock cycles for update\n", ((double) (end - start)));
			}	*/
		}

	}


	// change alpha beta osv osv

	for (int i = (*k); i>=2; --i)
	{
		alpha[i] = alpha[i-1];
	}
	alpha[0] = amount;
	alpha[1] = amount2;

	for (int i = ((*k)-1); i>=1; --i)
	{
		beta[i] = beta[i-1];
		gamma[i] = gamma[i-1];
	}

	beta[0] = 0;
	gamma[0]= 0;


	(*k)++;

	return;

}

// Reverse bottom split
void sparse_reverse_bs(dictionary * old_log_subperms,  dictionary * new_log_subperms,
						int amount, int amount2,
						int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m,
						int * k){

	int l=0;
	
	int r=0;

	double log_upper1=0;
	double log_upper2=0;
	double log_upper3=0;
	double log_lower1=0;
	double log_lower2=0;
	double log_lower3a=0;
	double log_lower3b=0;
	clock_t start, end;
	double value = 0;
	double new_value = 0;
	//double log_lower = 0;
	//
	pair new_pair = {.x = 0,.y=0};

	nullset_dictionary(new_log_subperms);

	for (int z = 0; z < (*old_log_subperms).used_len; ++z)
	{
		pair old_pair = (*old_log_subperms).array[z];

		r = old_pair.x;
		l = old_pair.y;

	
		value = (*old_log_subperms).value_array[z];

		// check bounds
		// 0<=r<=amount
		// max(l, m+amount + alpha[0] + alpha[k-1] - n -r)<= s <= min(alpha[k-1], l+amount)
		// r<=l<=r + amount2
		


		new_pair.x = r;

		for (int s = MAX(l, (*m)+amount + alpha[0] + alpha[(*k)-1] - n -r); s <= MIN(alpha[(*k)-1], l+amount); ++s)
		{
			// check bounds

			/*if(z){
				start = clock();
			}*/
			log_upper1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - ((*m) - r - l)];
			log_lower1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - (*m + amount) + r +s];
			log_upper2 = log_factorials[alpha[(*k)-1] -l  ];
			log_lower2 = log_factorials[alpha[(*k) -1] -s];
			log_upper3 = log_factorials[ amount ];
			log_lower3a = log_factorials[ amount - s + l ];
			log_lower3b = log_factorials[ s - l ];

			new_value = log_upper1 + log_upper2 + log_upper3 - log_lower1 - log_lower2 -
						   log_lower3a - log_lower3b + value;
			/*if(z<10){
				end = clock();
				Rprintf("used %f clock cycles for computing loglowerupperosv\n", ((double) (end - start)));
			}*/	
		
 			new_pair.y = s;
			/*if(z<10){
				start = clock();
			}*/
 			update_dict(new_pair, new_value, new_log_subperms);
 			/*if(z<10){
				end = clock();
				Rprintf("used %f clock cycles for update\n", ((double) (end - start)));
			}*/
		}
	}


	// change alpha beta osv osv
	
	gamma[0] = amount;
	*m = (*m) + amount;

	return;

}

// Reverse right merge
void sparse_reverse_rm(dictionary * old_log_subperms,  dictionary * new_log_subperms,
						int amount, int amount2,
						int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m,
						int * k){


	int r = 0;
	int l = 0;

	double log_upper1=0;
	double log_upper2=0;
	double log_upper3=0;
	double log_lower1=0;
	double log_lower2=0;
	double log_lower3=0;
	double log_upper4=0;
	double log_lower4a=0;
	double log_lower4b=0;


	double value = 0;
	double new_value = 0;
	//double log_lower = 0;
	//
	pair new_pair = {.x = 0,.y=0};

	nullset_dictionary(new_log_subperms);

	for (int z = 0; z < (*old_log_subperms).used_len; ++z)
	{
		pair old_pair = (*old_log_subperms).array[z];

		r = old_pair.x;
		l = old_pair.y;

		value = (*old_log_subperms).value_array[z];
		// check bounds
		// 0<=r<=amount
		// 0<= s <= alpha[(*k)-1]
		// r<=l<=r + amount2
		
		new_pair.x = r;

		for (int s = MAX(0, l-amount); s <= MIN(l, amount2); ++s)
		{

			log_upper1 = log_factorials[amount];
			log_upper2 = log_factorials[amount2];
			log_upper3 = log_factorials[amount + amount2];
			log_lower1 = log_factorials[amount - l +s];
			log_lower2 = log_factorials[amount2 - s];
			log_lower3 = log_factorials[amount +amount2 -l];
			log_upper4 = log_factorials[l];
			log_lower4a = log_factorials[l-s];
			log_lower4b = log_factorials[s];
			



			new_value = log_upper1 + log_upper2 - log_upper3 + log_upper4
								- log_lower1 - log_lower2 + log_lower3
								- log_lower4a -log_lower4b + value;
 			new_pair.y = s;

 			update_dict(new_pair, new_value, new_log_subperms);
		}
	}
	

	// change alpha beta osv osv
	
	alpha[(*k)-1] = amount;
	alpha[(*k)] = amount2;

	beta[(*k)-1] = 0;
	gamma[(*k)-1] = 0;

	
	(*k)++;

	return;

}

// Reverse top split
void sparse_reverse_ts(dictionary * old_log_subperms,  dictionary * new_log_subperms,
						int amount, int amount2,
						int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m,
						int * k){


	int l = 0;
	int s = 0;

	double log_upper1=0;
	double log_upper2=0;
	double log_upper3=0;
	double log_lower1=0;
	double log_lower2=0;
	double log_lower3a=0;
	double log_lower3b=0;


	double value = 0;
	double new_value = 0;
	//double log_lower = 0;
	//
	pair new_pair = {.x = 0,.y=0};

	nullset_dictionary(new_log_subperms);

	for (int z = 0; z < (*old_log_subperms).used_len; ++z)
	{
		pair old_pair = (*old_log_subperms).array[z];

		l = old_pair.x;
		s = old_pair.y;

		value = (*old_log_subperms).value_array[z];
		// check bounds
		// 0<=r<=amount
		// 0<= s <= alpha[(*k)-1]
		// r<=l<=r + amount2
		
		new_pair.y = s;

		for (int r = MAX(l,(*m)+amount +alpha[0] + alpha[(*k)-1] - n -s); r <= MIN(l + amount,alpha[0]); ++r)
		{

			
			log_upper1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - ((*m)   - l - s )];
			//log_lower1 = log_factorials[amount + l -s]
			log_lower1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - (*m + amount) + r +s];
			log_upper2 = log_factorials[alpha[0] - l];
			log_lower2 = log_factorials[alpha[0] - r];
			log_upper3 = log_factorials[ amount ];
			log_lower3a = log_factorials[ amount - r + l ];
			log_lower3b = log_factorials[ r-l];
			



			new_value = log_upper1 + log_upper2 + log_upper3 - log_lower1 - log_lower2 -
							   log_lower3a - log_lower3b + value;
 			new_pair.x = r;

 			update_dict(new_pair, new_value, new_log_subperms);
		}
	}
	

	// change alpha beta osv osv
	beta[(*k) - 2] = amount;
	*m = (*m) + amount;

	return;

}


void sparse_get_reduced_log_subperms(dictionary * new_log_subperms, int * alpha, int * beta, int * gamma, 
				double * log_factorials, int n, int * m, int * k){
	// case 1': 
	double new_value = 0;
	pair new_pair = {.x =0,.y=0};

	nullset_dictionary(new_log_subperms);


	if((*k)==2 && beta[0]==0 && gamma[0]==0){

		for (int r = 0; r <= (*m); ++r)
		{
			
			
			new_value= log_factorials[(*m)]- log_factorials[r] 
						- log_factorials[(*m)-r] + log_factorials[alpha[0]] - log_factorials[alpha[0]-r] 
						+ log_factorials[alpha[1]] - log_factorials[alpha[1]-((*m)-r)];
			
			new_pair.x = r;
			new_pair.y = (*m)-r;
			add_to_dictionary(new_log_subperms, new_pair, new_value);
			
		}
	}

	// case 2: 

	else if((*k)==2 && beta[0]==0 && gamma[0]==(*m)){
		new_pair.x = 0;
		new_pair.y = *m;

		new_value = log_factorials[alpha[1]] - log_factorials[alpha[1]-(*m)];
		add_to_dictionary(new_log_subperms, new_pair, new_value);

	}

	// case 3:
	else if((*k)==2 && beta[0]==(*m) && gamma[0]==0){
		new_pair.x = *m;
		new_pair.y = 0;

		new_value = log_factorials[alpha[0]] - log_factorials[alpha[0]-(*m)];
		add_to_dictionary(new_log_subperms, new_pair, new_value);

	}

	// case 4:
	else if((*k)==3 && beta[0]==0 && beta[1]==(*m) && gamma[0] == (*m) && gamma[1] == 0){
		new_pair.x = 0;
		new_pair.y = 0;
		new_value = log_factorials[alpha[1]] - log_factorials[alpha[1]-(*m)];
		add_to_dictionary(new_log_subperms, new_pair, new_value);
	}
	else{
		//Rprintf("Error in get_reduced_log_subperms! None of the four cases match the given arguments");
	}
}

dictionary* sparse_reverse_reduction(dictionary * old_log_subperms, dictionary * new_log_subperms, int * alpha, 
					   int * beta, int * gamma, int *m, int n, int *k, int * history,
			           int * amount_history, int * history_len, double * log_factorials){
	// INPUT OLD LOG SUBPERMS = the result of the reduced case
	
	// returns pointer to the new log subperms

	dictionary * tmp =NULL;

	int amount = 0;
	int amount2 = 0;

	for (int i =((*history_len)-1); i>=0; --i)
	{
		//Rprintf("i = %d\n",i);
		//Rprintf("old_subperms_len = %d\n", (*old_log_subperms).used_len);
		//Rprintf("old_subperms_maxsize = %d\n", (*old_log_subperms).table_size);
		//Rprintf("fraction =  %f\n",  ((double )(*old_log_subperms).used_len) / ((double) (n+1)*(n+1)));
		if(history[i]==0){
			// reverse top trim
			//Rprintf("Reverse top trim!\n");
			amount = amount_history[2*i];
			sparse_reverse_tt(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				k);
		}
		else if(history[i] ==1){
			// reverse bottom split
			//Rprintf("Reverse bottom split!\n");
			amount = amount_history[2*i];
			sparse_reverse_bs(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				 k);
		}

		else if(history[i]==2){
			// reverse left merge
			//Rprintf("Reverse left merge!\n");
			amount = amount_history[2*i];
			amount2 = amount_history[2*i+1];
			sparse_reverse_lm(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				 k);
		}

		else if(history[i]==3){
			// reverse bottom trim
			//Rprintf("Reverse bottom trim!\n");
			amount = amount_history[2*i];
			sparse_reverse_bt(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				k);
		}
		else if(history[i]==4){
			// reverse top split
			//Rprintf("Reverse top split!\n");
			amount = amount_history[2*i];
			sparse_reverse_ts(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				k);
		}
		else if(history[i]==5){
			// reverse right merge
			//Rprintf("Reverse right merge!\n");
			amount = amount_history[2*i];
			amount2 = amount_history[2*i+1];
			sparse_reverse_rm(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				k);
		}
		////Rprintf("old = %d\n", old_log_subperms);
		////Rprintf("new = %d\n", new_log_subperms);
		//Rprintf("historyindex = %d\n", i);
		//Rprintf("amount = %d\n", amount);
		////Rprintf("Result:\n");
		////Rprintf("Old:\n");
		//print_matrix(n+1,n+1,old_log_subperms);		
		////Rprintf("New:\n");
		//print_matrix(n+1,n+1,new_log_subperms);

		tmp  = old_log_subperms;
		old_log_subperms = new_log_subperms;
		new_log_subperms = tmp;
	}
	//print_matrix((n+1),(n+1), old_log_subperms);
	////Rprintf("\n");
	//print_matrix((n+1),(n+1), new_log_subperms);

	return old_log_subperms;
}


// 
// // Reverse bottom split
// void reverse_bs(double * old_log_subperms, double * new_log_subperms, int amount, int amount2,
// 				int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m, 
// 				double * temp_vec, int * k){
// 
// 	for (int i = 0; i < (n+1)*(n+1); ++i)
// 	{
// 		new_log_subperms[i]=-1;
// 	}
// 
// 	////Rprintf("HEIIII\n");
// 	double log_upper1=0;
// 	double log_upper2=0;
// 	double log_upper3=0;
// 	double log_lower1=0;
// 	double log_lower2=0;
// 	double log_lower3a=0;
// 	double log_lower3b=0;
// 
// 
// 	double maxval = -1;
// 	int lstart = 0;
// 
// 	////Rprintf("alphas = ");
// 	//print_int_vector(n,alpha);
// 	////Rprintf("n - alpha[0] - alpha[k-1] = %d\n",n - alpha[0] - alpha[(*k)-1]  );
// 	for (int r = 0; r <= alpha[0]; ++r)
// 	{
// 		for (int s = MAX(0, (*m + amount) +alpha[0] + alpha[(*k)-1]-n-r); s<= alpha[(*k)-1]; ++s)
// 		{
// 			maxval=-1;
// 			lstart = MAX(s-amount,0);
// 			for (int l = lstart; l <= s; ++l)
// 			{
// 
// 				////Rprintf("r = %d, s = %d, l = %d\n", r,s,l);
// 
// 				if(r == 1 && s == 1){
// 					////Rprintf("hit r=1, s=1, l=%d!\n", l);
// 					////Rprintf("old_log_subperms[r,l] = %f\n", old_log_subperms[cord_spec(r,l,(n+1))]);
// 
// 				}
// 
// 				if(old_log_subperms[cord_spec(r,l,(n+1))] < 0 ){
// 					////Rprintf("SKIPPED\n");
// 					temp_vec[l-lstart] = -1;
// 					continue;
// 				}
// 
// 
// 
// 
// 				log_upper1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - ((*m) - r - l)];
// 				//log_lower1 = log_factorials[amount + l -s]
// 				log_lower1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - (*m + amount) + r +s];
// 				log_upper2 = log_factorials[alpha[(*k)-1] -l  ];
// 				log_lower2 = log_factorials[alpha[(*k) -1] -s];
// 				log_upper3 = log_factorials[ amount ];
// 				log_lower3a = log_factorials[ amount - s + l ];
// 				log_lower3b = log_factorials[ s - l ];
// 
// 				temp_vec[l-lstart] = log_upper1 + log_upper2 + log_upper3 - log_lower1 - log_lower2 -
// 							   log_lower3a - log_lower3b + old_log_subperms[cord_spec(r,l,(n+1))];
// 
// 				if(r==1 && s==1){
// 					////Rprintf("log_upper1 = %f\n", log_upper1);
// 					////Rprintf("log_upper2 = %f\n", log_upper2);
// 					////Rprintf("log_upper3 = %f\n", log_upper3);
// 					////Rprintf("log_lower1 = %f\n", log_lower1);
// 					////Rprintf("log_lower2 = %f\n", log_lower2);
// 					////Rprintf("log_lower3a = %f\n", log_lower3a);
// 					////Rprintf("log_lower3b = %f\n", log_lower3b);
// 
// 				}
// 
// 				if(temp_vec[l-lstart] > maxval){
// 					maxval = temp_vec[l-lstart];
// 				}
// 
// 
// 
// 			}
// 
// 			// log sum exp her
// 
// 			new_log_subperms[cord_spec(r,s,(n+1))] = Clog_sum_exp(temp_vec, s-lstart+1, maxval);
// 			
// 		}
// 	}
// 
// 	// change alpha beta osv osv
// 	
// 	gamma[0] = amount;
// 	*m = (*m) + amount;
// 
// 	return;
// 	
// }
// 
// // Reverse top split
// void reverse_ts(double * old_log_subperms, double * new_log_subperms, int amount, int amount2,
// 				int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m, 
// 				double * temp_vec, int * k){
// 
// 	for (int i = 0; i < (n+1)*(n+1); ++i)
// 	{
// 		new_log_subperms[i]=-1;
// 	}
// 
// 	double log_upper1=0;
// 	double log_upper2=0;
// 	double log_upper3=0;
// 	double log_lower1=0;
// 	double log_lower2=0;
// 	double log_lower3a=0;
// 	double log_lower3b=0;
// 
// 
// 	double maxval = -1;
// 	int lstart = 0;
// 
// 	for (int r = 0; r <= alpha[0]; ++r)
// 	{
// 		for (int s = MAX((*m+amount+alpha[0]+alpha[(*k)-1] -n -r),0); s<= alpha[(*k)-1]; ++s)
// 		{
// 			maxval=-1;
// 		  lstart = MAX(0, r - amount);
// 			for (int l = lstart; l <= r; ++l)
// 			{
// 				if(old_log_subperms[cord_spec(l,s,(n+1))] < 0 ){
// 					temp_vec[l-lstart] = -1;
// 					continue;
// 				}
// 
// 
// 				log_upper1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - ((*m)   - l - s )];
// 				//log_lower1 = log_factorials[amount + l -s]
// 				log_lower1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - (*m + amount) + r +s];
// 				log_upper2 = log_factorials[alpha[0] - l];
// 				log_lower2 = log_factorials[alpha[0] - r];
// 				log_upper3 = log_factorials[ amount ];
// 				log_lower3a = log_factorials[ amount - r + l ];
// 				log_lower3b = log_factorials[ r-l];
// 
// 				temp_vec[l-lstart] = log_upper1 + log_upper2 + log_upper3 - log_lower1 - log_lower2 -
// 							   log_lower3a - log_lower3b + old_log_subperms[cord_spec(l,s,(n+1))];
// 
// 
// 				if(temp_vec[l-lstart] > maxval){
// 					maxval = temp_vec[l-lstart];
// 				}
// 
// 
// 
// 			}
// 
// 			// log sum exp her
// 
// 			new_log_subperms[cord_spec(r,s,(n+1))] = Clog_sum_exp(temp_vec, r - lstart+1, maxval);
// 			
// 		}
// 	}
// 
// 	// change alpha beta osv osv
// 	beta[(*k) - 2] = amount;
// 	*m = (*m) + amount;
// 
// 	return;
// 	
// }
// 
// 
// // Reverse left merge
// void reverse_lm(double * old_log_subperms, double * new_log_subperms, int amount, int amount2,
// 				int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m, 
// 				double * temp_vec, int * k){
// 
// 
// 	// amount = alpha[0]
// 	// amount2 = alpha[1]
// 
// 	for (int i = 0; i < (n+1)*(n+1); ++i)
// 	{
// 		new_log_subperms[i]=-1;
// 	}
// 
// 	double log_upper1=0;
// 	double log_upper2=0;
// 	double log_upper3=0;
// 	double log_lower1=0;
// 	double log_lower2=0;
// 	double log_lower3=0;
// 	double log_upper4=0;
// 	double log_lower4a=0;
// 	double log_lower4b=0;
// 
// 
// 	double maxval = -1;
// 
// 	log_upper1 = log_factorials[amount2];
// 	log_upper2 = log_factorials[amount];
// 	log_upper3 = log_factorials[amount + amount2];
// 
// for (int r = 0; r <= (amount); ++r)
// {
// 	log_lower2 = log_factorials[amount - r];
// 	for (int s = 0; s<= alpha[(*k)-1]; ++s)
// 	{
// 		maxval=-1;
// 		for (int l = r; l <= (r+amount2); ++l)
// 		{
// 			////Rprintf("r = %d, s = %d, l = %d\n",r,s,l);
// 			if(old_log_subperms[cord_spec(l,s,(n+1))] < 0 ){
// 				temp_vec[l-r] = -1;
// 				////Rprintf("skipping!\n");
// 				continue;
// 				}
// 
// 
// 				log_lower1 = log_factorials[amount2 - l +r];
// 				log_lower3 = log_factorials[ amount +amount2 -l];
// 				log_upper4 = log_factorials[l];
// 				log_lower4a = log_factorials[l-r];
// 				log_lower4b = log_factorials[r];
// 

// 
// 				temp_vec[l-r] = log_upper1 + log_upper2 - log_upper3 + log_upper4
// 								- log_lower1 - log_lower2 + log_lower3
// 								- log_lower4a -log_lower4b + old_log_subperms[cord_spec(l,s,(n+1))];
// 
// 
// 				if(temp_vec[l-r] > maxval){
// 					maxval = temp_vec[l-r];
// 				}
// 
// 
// 			}
// 
// 			// log sum exp her
// 
// 			new_log_subperms[cord_spec(r,s,(n+1))] = Clog_sum_exp(temp_vec, (amount2+1), maxval);
// 			
// 		}
// 	}
// 
	// // change alpha beta osv osv
	// 
	// for (int i = (*k); i>=2; --i)
	// {
	// 	alpha[i] = alpha[i-1];
	// }
	// alpha[0] = amount;
	// alpha[1] = amount2;
	// 
	// for (int i = ((*k)-1); i>=1; --i)
	// {
	// 	beta[i] = beta[i-1];
	// 	gamma[i] = gamma[i-1];
	// }
	// 
	// beta[0] = 0;
	// gamma[0]= 0;
	// 
	// 
	// 
	// (*k)++;
	// 
	// return;
// 	
// }
// 
// 
// 
// // Reverse right merge
// void reverse_rm(double * old_log_subperms, double * new_log_subperms, int amount, int amount2,
// 				int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m, 
// 				double * temp_vec, int * k){
// 
// 	for (int i = 0; i < (n+1)*(n+1); ++i)
// 	{
// 		new_log_subperms[i]=-1;
// 	}
// 
// 	// amount = alpha[0]
// 	// amount2 = alpha[1]
// 
// 
// 	double log_upper1=0;
// 	double log_upper2=0;
// 	double log_upper3=0;
// 	double log_lower1=0;
// 	double log_lower2=0;
// 	double log_lower3=0;
// 	double log_upper4=0;
// 	double log_lower4a=0;
// 	double log_lower4b=0;
// 
// 
// 	double maxval = -1;
// 
// 	log_upper1 = log_factorials[amount];
// 	log_upper2 = log_factorials[amount2];
// 	log_upper3 = log_factorials[amount + amount2];
// 
// 	for (int r = 0; r <= alpha[0]; ++r)
// 	{
// 		log_lower2 = log_factorials[amount - r];
// 		for (int s = 0; s<= (amount2); ++s)
// 		{
// 			maxval=-1;
// 			for (int l = s; l <= (amount+s); ++l)
// 			{
// 
// 				if(old_log_subperms[cord_spec(r,l,(n+1))] < 0 ){
// 					temp_vec[l-s] = -1;
// 					continue;
// 				}
// 
// 
// 				log_lower1 = log_factorials[amount - l +s];
// 				log_lower2 = log_factorials[amount2 - s];
// 				log_lower3 = log_factorials[amount +amount2 -l];
// 				log_upper4 = log_factorials[l];
// 				log_lower4a = log_factorials[l-s];
// 				log_lower4b = log_factorials[s];
// 
// 			
// 
// 				temp_vec[l-s] = log_upper1 + log_upper2 - log_upper3 + log_upper4
// 								- log_lower1 - log_lower2 + log_lower3
// 								- log_lower4a -log_lower4b + old_log_subperms[cord_spec(r,l,(n+1))];
// 
// 
// 				if(temp_vec[l-s] > maxval){
// 					maxval = temp_vec[l-s];
// 				}
// 
// 
// 			}
// 
// 			// log sum exp her
// 
// 			new_log_subperms[cord_spec(r,s,(n+1))] = Clog_sum_exp(temp_vec, (amount+1), maxval);
// 			
// 		}
// 	}
// 
// 	// change alpha beta osv osv
// 	
// 	alpha[(*k)-1] = amount;
// 	alpha[(*k)] = amount2;
// 
// 	beta[(*k)-1] = 0;
// 	gamma[(*k)-1] = 0;
// 
// 	
// 	(*k)++;
// 
// 	return;
// }
// 
// 
// void get_reduced_log_subperms(double * new_log_subperms, int * alpha, int * beta, int * gamma, 
// 				double * log_factorials, int n, int * m, int * k){
// 	// case 1': 
// 	if((*k)==2 && beta[0]==0 && gamma[0]==0){
// 
// 		for (int r = 0; r < (n+1); ++r)
// 		{
// 			for (int s = 0; s < (n+1); ++s)
// 			{
// 				if(r+s==(*m)){
// 					new_log_subperms[cord_spec(r,s,(n+1))] = log_factorials[(*m)]- log_factorials[r] 
// 					- log_factorials[(*m)-r] + log_factorials[alpha[0]] - log_factorials[alpha[0]-r] 
// 					+ log_factorials[alpha[1]] - log_factorials[alpha[1]-s];
// 				}
// 				else{
// 					new_log_subperms[cord_spec(r,s,(n+1))] = -1;
// 				}
// 			}
// 		}
// 	}
// 
// 	// case 2: 
// 
// 	else if((*k)==2 && beta[0]==0 && gamma[0]==(*m)){
// 		for (int r = 0; r < (n+1); ++r)
// 		{
// 			for (int s = 0; s < (n+1); ++s)
// 			{
// 				
// 				new_log_subperms[cord_spec(r,s,(n+1))] = -1;
// 				
// 			}
// 		}
// 		new_log_subperms[cord_spec(0,(*m), (n+1))] = log_factorials[alpha[1]] - log_factorials[alpha[1]-(*m)];
// 
// 	}
// 
// 	// case 3:
// 	else if((*k)==2 && beta[0]==(*m) && gamma[0]==0){
// 		for (int r = 0; r < (n+1); ++r)
// 		{
// 			for (int s = 0; s < (n+1); ++s)
// 			{
// 				
// 				new_log_subperms[cord_spec(r,s,(n+1))] = -1;
// 				
// 			}
// 		}
// 		new_log_subperms[cord_spec((*m),0,(n+1))] = log_factorials[alpha[0]] - log_factorials[alpha[0]-(*m)];
// 	}
// 
// 	// case 4:
// 	else if((*k)==3 && beta[0]==0 && beta[1]==(*m) && gamma[0] == (*m) && gamma[1] == 0){
// 		for (int r = 0; r < (n+1); ++r)
// 		{
// 			for (int s = 0; s < (n+1); ++s)
// 			{
// 				
// 				new_log_subperms[cord_spec(r,s,(n+1))] = -1;
// 				
// 			}
// 		}
// 		new_log_subperms[cord_spec(0,0,(n+1))] = log_factorials[alpha[1]] - log_factorials[alpha[1]-(*m)];
// 	}
// 	else{
// 		//Rprintf("Error in get_reduced_log_subperms! None of the four cases match the given arguments");
// 	}
// }
// 
// char check_if_reduced(int * alpha, int * beta, int * gamma, 
// 					  int * m, int * k){
// 	// case 1': 
// 	char reduced = 0;
// 
// 	if((*k)==2 && beta[0]==0 && gamma[0]==0){
// 		reduced =1;
// 	}
// 
// 	// case 2: 
// 
// 	else if((*k)==2 && beta[0]==0 && gamma[0]==(*m)){
// 		reduced = 1;
// 
// 	}
// 
// 	// case 3:
// 	else if((*k)==2 && beta[0]==(*m) && gamma[0]==0){
// 		reduced = 1;
// 	}
// 
// 
// 	// case 4:
// 	else if((*k)==3 && beta[0]==0 && beta[1]==(*m) && gamma[0] == (*m) && gamma[1] == 0){
// 		reduced = 1;
// 	}
// 	return reduced;
// }
// 
// 
// void reduction(int * alpha, int * beta, int * gamma, int *m, int n, int *k, int * history,
// 			   int * amount_history, int * history_len){
// 
// 	// legg inn kode i tilfelle history eller amounthistry == NULL
// 	
// 	// assume history and amount_history are of size 3n and 6n, respectively
// 	// amount_history assumed 3n \times 2 matrix, column wise stored
// 
// 	// history codes:
// 	// 0 = top trim
// 	// 1 = bottom split
// 	// 2 = left merge
// 	// 3 = bottom trim
// 	// 4 = top split
// 	// 5 = right merge
// 
// 	memset(history, 0, sizeof(int)*3*n );
// 	memset(amount_history, 0, sizeof(int)*6*n );
// 	*history_len =0;
// 
// 
// 	char reduced = check_if_reduced(alpha, beta, gamma, 
// 					  m, k);
// 
// 	int historycounter = 0;
// 
// 	int cc = 0;
// 
// 
// 	while(1){
// 		// if top trim possible
// 		if((*m)>beta[0] && beta[0] > 0){
// 			history[historycounter] = 0;
// 			amount_history[2*historycounter] = beta[0];
// 			(*m) = (*m)-beta[0];
// 			beta[0]=0;
// 			historycounter++;
// 			//Rprintf("TOP TRIM ->\n");
// 			//print_int_vector(n, beta);
// 			//Rprintf("=====\n");
// 			//Rprintf("alpha:\n");
// 			print_int_vector(n, alpha);
// 			//Rprintf("beta:\n");
// 			print_int_vector(n, beta);
// 			//Rprintf("gamma:\n");
// 			print_int_vector(n, gamma);
// 			//Rprintf("m:\n");
// 			print_int_vector(1, m);
// 			//Rprintf("k:\n");
// 			print_int_vector(1, k);
// 			//Rprintf("=====\n");
// 			if(check_if_reduced(alpha, beta, gamma, 
// 					  m, k)){break;}
// 
// 		}
// 		// if bottom split possible
// 		else if((*m)>gamma[0] && gamma[0] > 0){
// 			history[historycounter] = 1;
// 			amount_history[2*historycounter] = gamma[0];
// 			(*m) = (*m) - gamma[0];
// 			gamma[0] = 0;
// 			historycounter++;
// 			//Rprintf("BOTTOM SPLIT ->\n");
// 			//print_int_vector(n, beta);
// 			//Rprintf("=====\n");
// 			//Rprintf("alpha:\n");
// 			print_int_vector(n, alpha);
// 			//Rprintf("beta:\n");
// 			print_int_vector(n, beta);
// 			//Rprintf("gamma:\n");
// 			print_int_vector(n, gamma);
// 			//Rprintf("m:\n");
// 			print_int_vector(1, m);
// 			//Rprintf("k:\n");
// 			print_int_vector(1, k);
// 			//Rprintf("=====\n");
// 			if(check_if_reduced(alpha, beta, gamma, 
// 					  m, k)){break;}
// 		}
// 
// 		// if bottom left merge
// 		else if(beta[0]==0 && gamma[0]==0){
// 			history[historycounter] = 2;
// 			amount_history[2*historycounter] = alpha[0];
// 			amount_history[2*historycounter+1] = alpha[1];
// 
// 			//shift alpha:
// 			alpha[0] = alpha[0]+ alpha[1];
// 			for (int i = 1; i < (*k-1); ++i)
// 			{
// 				alpha[i] = alpha[i+1];
// 			}
// 
// 			for (int i = 0; i < (*k-2); ++i)
// 			{
// 				beta[i] = beta[i+1];
// 				gamma[i] = gamma[i+1];
// 			}
// 			(*k)--;
// 			historycounter++;
// 			//Rprintf("LEFT MERGE ->\n");
// 			//print_int_vector(n, beta);
// 			//Rprintf("=====\n");
// 			//Rprintf("alpha:\n");
// 			print_int_vector(n, alpha);
// 			//Rprintf("beta:\n");
// 			print_int_vector(n, beta);
// 			//Rprintf("gamma:\n");
// 			print_int_vector(n, gamma);
// 			//Rprintf("m:\n");
// 			print_int_vector(1, m);
// 			//Rprintf("k:\n");
// 			print_int_vector(1, k);
// 			//Rprintf("=====\n");
// 			if(check_if_reduced(alpha, beta, gamma, 
// 					  m, k)){break;}
// 		}		
// 		// reduced = check_if_reduced(alpha, beta, gamma, 
// 		// 			  m, k);
// 		/*if(reduced){
// 			break;
// 		}*/
// 		// if bottom bottom trim
// 		else if((*m) > gamma[(*k)-2] && gamma[(*k)-2]>0){
// 			history[historycounter] = 3;
// 			amount_history[2*historycounter] = gamma[(*k)-2];
// 
// 			(*m) = (*m) - gamma[(*k)-2];
// 			gamma[(*k)-2]=0;
// 			historycounter++;
// 			//Rprintf("BOTTOM TRIM ->\n");
// 			//print_int_vector(n, beta);
// 			//Rprintf("=====\n");
// 			//Rprintf("alpha:\n");
// 			print_int_vector(n, alpha);
// 			//Rprintf("beta:\n");
// 			print_int_vector(n, beta);
// 			//Rprintf("gamma:\n");
// 			print_int_vector(n, gamma);
// 			//Rprintf("m:\n");
// 			print_int_vector(1, m);
// 			//Rprintf("k:\n");
// 			print_int_vector(1, k);
// 			//Rprintf("=====\n");
// 			if(check_if_reduced(alpha, beta, gamma, 
// 					  m, k)){break;}
// 		}
// 
// 		// if top split possible
// 		else if( (*m)> beta[(*k)-2] && beta[(*k)-2]>0 ){
// 			history[historycounter] = 4;
// 			amount_history[2*historycounter] = beta[(*k)-2];
// 			(*m) = (*m) - beta[(*k)-2];
// 			beta[(*k)-2]=0;
// 			historycounter++;
// 			//Rprintf("TOP SPLIT ->\n");
// 			//print_int_vector(n, beta);
// 			//Rprintf("=====\n");
// 			//Rprintf("alpha:\n");
// 			print_int_vector(n, alpha);
// 			//Rprintf("beta:\n");
// 			print_int_vector(n, beta);
// 			//Rprintf("gamma:\n");
// 			print_int_vector(n, gamma);
// 			//Rprintf("m:\n");
// 			print_int_vector(1, m);
// 			//Rprintf("k:\n");
// 			print_int_vector(1, k);
// 			//Rprintf("=====\n");
// 			if(check_if_reduced(alpha, beta, gamma, 
// 					  m, k)){break;}
// 		}
// 
// 		// if right merge possible
// 		else if(beta[(*k)-2]==0 && gamma[(*k)-2]==0 ){
// 			history[historycounter] = 5;
// 			amount_history[2*historycounter] = alpha[(*k)-2];
// 			amount_history[2*historycounter+1] = alpha[(*k)-1];
// 
// 
// 			alpha[(*k)-2] = alpha[(*k)-2] + alpha[(*k)-1];
// 			(*k)--;
// 
// 			historycounter++;
// 			//Rprintf("RIGHT MERGE ->\n");
// 			//print_int_vector(n, beta);
// 			//Rprintf("=====\n");
// 			//Rprintf("alpha:\n");
// 			print_int_vector(n, alpha);
// 			//Rprintf("beta:\n");
// 			print_int_vector(n, beta);
// 			//Rprintf("gamma:\n");
// 			print_int_vector(n, gamma);
// 			//Rprintf("m:\n");
// 			print_int_vector(1, m);
// 			//Rprintf("k:\n");
// 			print_int_vector(1, k);
// 			//Rprintf("=====\n");
// 			if(check_if_reduced(alpha, beta, gamma, 
// 					  m, k)){break;}
// 		}		
// 
// 		cc++;
// 
// 
// 	}
// 
// 	//Rprintf("Terminated at iter %d\n", cc);
// 	
// 	*history_len = historycounter;
// }
// 
// double * reverse_reduction(double * old_log_subperms, double * new_log_subperms, int * alpha, 
// 					   int * beta, int * gamma, int *m, int n, int *k, int * history,
// 			           int * amount_history, int * history_len, double * log_factorials, 
// 			           double * temp_vec){
// 	// INPUT OLD LOG SUBPERMS = the result of the reduced case
// 	
// 	// returns pointer to the new log subperms
// 
// 	double * tmp =NULL;
// 
// 	int amount = 0;
// 	int amount2 = 0;
// 
// 	for (int i =((*history_len)-1); i>=0; --i)
// 	{
// 		
// 		if(history[i]==0){
// 			// reverse top trim
// 			//Rprintf("Reverse top trim!\n");
// 			amount = amount_history[2*i];
// 			reverse_tt(old_log_subperms, new_log_subperms, amount, amount2,
// 				alpha,  beta, gamma,log_factorials, n, m, 
// 				k);
// 		}
// 		else if(history[i] ==1){
// 			// reverse bottom split
// 			//Rprintf("Reverse bottom split!\n");
// 			amount = amount_history[2*i];
// 			reverse_bs(old_log_subperms, new_log_subperms, amount, amount2,
// 				alpha,  beta, gamma,log_factorials, n, m, 
// 				temp_vec, k);
// 		}
// 
// 		else if(history[i]==2){
// 			// reverse left merge
// 			//Rprintf("Reverse left merge!\n");
// 			amount = amount_history[2*i];
// 			amount2 = amount_history[2*i+1];
// 			reverse_lm(old_log_subperms, new_log_subperms, amount, amount2,
// 				alpha,  beta, gamma,log_factorials, n, m, 
// 				temp_vec, k);
// 		}
// 
// 		else if(history[i]==3){
// 			// reverse bottom trim
// 			//Rprintf("Reverse bottom trim!\n");
// 			amount = amount_history[2*i];
// 			reverse_bt(old_log_subperms, new_log_subperms, amount, amount2,
// 				alpha,  beta, gamma,log_factorials, n, m, 
// 				k);
// 		}
// 		else if(history[i]==4){
// 			// reverse top split
// 			//Rprintf("Reverse top split!\n");
// 			amount = amount_history[2*i];
// 			reverse_ts(old_log_subperms, new_log_subperms, amount, amount2,
// 				alpha,  beta, gamma,log_factorials, n, m, 
// 				temp_vec,k);
// 		}
// 		else if(history[i]==5){
// 			// reverse right merge
// 			//Rprintf("Reverse right merge!\n");
// 			amount = amount_history[2*i];
// 			amount2 = amount_history[2*i+1];
// 			reverse_rm(old_log_subperms, new_log_subperms, amount, amount2,
// 				alpha,  beta, gamma,log_factorials, n, m, 
// 				temp_vec,k);
// 		}
// 		////Rprintf("old = %d\n", old_log_subperms);
// 		////Rprintf("new = %d\n", new_log_subperms);
// 		//Rprintf("historyindex = %d\n", i);
// 		//Rprintf("amount = %d\n", amount);
// 		////Rprintf("Result:\n");
// 		////Rprintf("Old:\n");
// 		//print_matrix(n+1,n+1,old_log_subperms);		
// 		////Rprintf("New:\n");
// 		//print_matrix(n+1,n+1,new_log_subperms);
// 
// 		tmp  = old_log_subperms;
// 		old_log_subperms = new_log_subperms;
// 		new_log_subperms = tmp;
// 	}
// 	//print_matrix((n+1),(n+1), old_log_subperms);
// 	////Rprintf("\n");
// 	//print_matrix((n+1),(n+1), new_log_subperms);
// 
// 	return old_log_subperms;
// }
// 


