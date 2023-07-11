#include "header.h"


size_t hash(pair p, int table_size) {
    // A simple hash function that combines the x and y values
    int t[2];
    t[0] = p.x;
    t[1] = p.y;
    XXH64_hash_t hash = XXH64(t, sizeof(int)*2, 1337);
    //Rprintf("%lx\n", hash);
/*    Rprintf("%ld\n", (int)hash);
    Rprintf("%ld\n", abs((int)hash));
    Rprintf("%d\n", sizeof(size_t));*/
    //Rprintf("hash(%d,%d): %d\n",p.x, p.y, abs(((int) hash) % table_size));
    return( hash % table_size);
    //return (4294967279 * p.x + 4294967291 * p.y) % table_size;
    //return(1001*p.x + p.y);
    
    //return ((p.x ^ p.y) *4294967279 ) % table_size;

}



double lookup_logperm(dictionary * dict, pair p){
	size_t h = hash(p, (*dict).table_size);
	//int index = -1;
	while ((*dict).table[h] >= 0) {
        if ((*dict).array[(*dict).table[h]].x == p.x && (*dict).array[(*dict).table[h]].y == p.y) {
            //index = dict.table[h];
            return (*dict).value_array[(*dict).table[h]];
            //break;
        }
        h = (h + 1) % (*dict).table_size;
    }

    return -1;
}


void add_to_dictionary(dictionary * dict, pair p, double val){
	// antar at vi VET at (key,value) ikke er i dict


	if ((RATIO*(*dict).used_len +1)> (*dict).table_size)
	{	
		
		expand_dictionary(dict);
	}

	size_t used_len = (*dict).used_len;

	//(*dict).array[used_len] = p;
	size_t h = hash(p, (*dict).table_size);
    while ((*dict).table[h] >= 0) {
            h = (h + 1) % (*dict).table_size;  // Linear probing to find the next empty slot
    }
    (*dict).table[h] = used_len;
	(*dict).value_array[used_len] = val;
	(*dict).array[used_len] = p;
	(*dict).used_len++;

}

void update_dict(pair p, double value, dictionary * dict){

	// check if (r,s) is already key in dict:
	
	size_t h = hash(p, (*dict).table_size);
	int collcount = 0;
	//int index = -1;
	int index = -1;
	while ((*dict).table[h] >= 0) {
		collcount++;
        if ((*dict).array[(*dict).table[h]].x == p.x && (*dict).array[(*dict).table[h]].y == p.y) {
            //index = dict.table[h];
            index = (*dict).table[h];
            break;
        }
        h = (h + 1) % (*dict).table_size;
        //collcount++;
    }

/*    
    if(collcount>10){
    	Rprintf("collcount = %d\n",collcount);
    }
*/    
    if(index ==-1 || index >=(*dict).used_len){
    	//if(collcount>10){
		//    Rprintf("[%d,%d]-> h = %ld : collcount = %d, table_size = %ld, used_len = %ld\n", p.x, p.y,h,collcount, (*dict).table_size, (*dict).used_len);
    	//}
    	// NEW VALUE!!!
    	if (RATIO*((*dict).used_len +1)> (*dict).table_size)
		{	
			//Rprintf("EXPANDING DICTIONARY: used = %ld, table_size = %ld\n",(*dict).used_len +1,(*dict).table_size);
			expand_dictionary(dict);
		}

		size_t used_len = (*dict).used_len;
	    (*dict).table[h] = used_len;
		(*dict).value_array[used_len] = value;
		(*dict).array[used_len] = p;
		(*dict).used_len++;
		return;
    }

    double existing_val = (*dict).value_array[index];
    //Rprintf("FOUND IN DICT! existing_val = %f\n", existing_val);
    //Rprintf("new input value = %f\n", value);

    if(existing_val>value){
    	//Rprintf("existing_val > value\n");
    	(*dict).value_array[index] += log(1.0 + exp(value - existing_val));
    }else{
    	(*dict).value_array[index] = value + log(1.0+ exp(existing_val - value));
    }



}



dictionary * init_dictionary(size_t init_size){
	dictionary * dict = (dictionary *) R_Calloc(1, dictionary);
	(*dict).array = (pair *) R_Calloc(init_size,  pair);
	(*dict).used_len = 0;
	(*dict).table_size = init_size;
	(*dict).value_array = (double*) R_Calloc(init_size, double);
	(*dict).table = (int*) R_Calloc(init_size, int);

	for (size_t i = 0; i < (*dict).table_size; i++) {
        (*dict).table[i] = -1;  // Initialize all entries to -1 (empty)
    }
	return dict;
}

void free_dictionary(dictionary * dict){
	R_Free((*dict).array);
	R_Free((*dict).value_array);
	R_Free(dict);
}

void expand_dictionary(dictionary * dict){
	//Rprintf("EXPANDING DICTIONARY. old_size = %ld, new = %ld\n", (*dict).table_size, 2*(*dict).table_size);
	pair * newarray = R_Realloc((*dict).array, (size_t) 2*(*dict).table_size, pair);
	if(newarray==NULL){
		newarray= (pair *) R_Calloc( 2 * ((*dict).table_size),  pair);
		memcpy(newarray, (*dict).array, sizeof(pair)*(*dict).used_len);
		R_Free((*dict).array);
	}
	(*dict).array = newarray;
	
	//
	//double * newvalue_array = (double *) R_Calloc( 2 * ((*dict).table_size),  double);
	double * newvalue_array = (double *) R_Realloc((*dict).value_array, (size_t) 2 * ((*dict).table_size),  double);
	if(newvalue_array==NULL){
		newvalue_array= (double *) R_Calloc( 2 * ((*dict).table_size),  double);
		memcpy(newvalue_array, (*dict).value_array, sizeof(double)*(*dict).used_len);
		R_Free((*dict).value_array);
	}
		
	(*dict).value_array = newvalue_array;
	

	int * newtable = (int *) R_Realloc((*dict).table, (size_t) 2 * ((*dict).table_size),  int);
	if(newtable==NULL){
		newtable= (int *) R_Calloc( 2 * ((*dict).table_size),  int);
		memcpy(newtable, (*dict).table, sizeof(int)*(*dict).table_size);
		R_Free((*dict).table);
		
	}
	(*dict).table = newtable;
	
	(*dict).table_size = 2*(*dict).table_size;


	for (int i = 0; i < (*dict).table_size; i++) {
        (*dict).table[i] = -1;  // Initialize all entries to -1 (empty)
    }
    for (int i = 0; i < (*dict).used_len; i++) {
        size_t h = hash((*dict).array[i], (*dict).table_size);
        while ((*dict).table[h] >= 0) {
            h = (h + 1) % (*dict).table_size;  // Linear probing to find the next empty slot
        }
        (*dict).table[h] = i;
    }

}

void nullset_dictionary(dictionary * dict){
	(*dict).used_len = 0;
	for (int i = 0; i < (*dict).table_size; i++) {
        (*dict).table[i] = -1;  // Initialize all entries to -1 (empty)
    }
}




void print_sparse_matrix(dictionary * matrix) {

    // Print out the matrix with aligned columns
    for (int i = 0; i < (*matrix).used_len; i++) {
    	Rprintf("(%d, %d) = %f\n", (*matrix).array[i].x, (*matrix).array[i].y, (*matrix).value_array[i]);

    }
}


void print_matrix(int rows, int cols, double * matrix) {
	return;

    // Print out the matrix with aligned columns
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            Rprintf("%5f ", matrix[cord_spec(i,j,rows)]);
        }
        //Rprintf("\n");
    }
}

void print_int_vector(int len,  int * vec) {
	/*return;*/

    // Print out the matrix with aligned columns
    Rprintf("(");
    for (int i = 0; i < (len); i++) {
        Rprintf("%d ", vec[i]);
    }
    Rprintf(")\n");
}
void print_float_vector(int len,  double * vec) {
	/*return;*/

    // Print out the matrix with aligned columns
    Rprintf("(");
    for (int i = 0; i < (len); i++) {
        Rprintf("%f ", vec[i]);
    }
    Rprintf(")\n");
}



double Clog_sum_exp(double * array, int len, double maxval){

	// ignore neg values

	if(maxval<0){
		return -1;
	}

	double exp_result = 0;



	for (int i = 0; i < len; ++i)
	{
		if(array[i]<0){
			continue;
		}

		exp_result += exp(array[i] - maxval);
	}

	////Rprintf("res = %f\n", (maxval + log(exp_result)));
	return (maxval + log(exp_result));
}

double Csparse_log_sum_exp(dictionary * dict){

	// ignore neg values

	double maxval = -1;

	for (int z = 0; z < (*dict).used_len; ++z)
	{
		if((*dict).value_array[z]> maxval){
			 maxval = (*dict).value_array[z];
		}
	}

	if(maxval<0){
		return -1;
	}

	double exp_result = 0;



	for (int z = 0; z < (*dict).used_len; ++z)
	{

		exp_result += exp((*dict).value_array[z] - maxval);
	}

	////Rprintf("res = %f\n", (maxval + log(exp_result)));
	return (maxval + log(exp_result));
}

SEXP C_wrapper_log_sum_exp(SEXP xSEXP, SEXP lenSEXP){

	// ignore neg values

	PROTECT(xSEXP);
	PROTECT(lenSEXP);

	int len = *(INTEGER(lenSEXP));
	double *x = REAL(xSEXP);

	SEXP retSEXP = PROTECT(allocVector(REALSXP, 1));
	double * ret = REAL(retSEXP);

	double maxval = -1;
	for (int i = 0; i < len; ++i)
	{
		if(x[i]> maxval){
			maxval = x[i];
		}
	}
	////Rprintf("maxval = %f\n", maxval);
	*ret = Clog_sum_exp(x, len, maxval);

	UNPROTECT(3);
	return retSEXP;


}
