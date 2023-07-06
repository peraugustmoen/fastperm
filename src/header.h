#define cord(r,c) ((r) + (DEPTH)*(c))
#define cord_spec(r,c, D) ((r) + (D)*(c))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
/*#include <R_ext/RS.h>	// RK addition
#include <R_ext/Lapack.h> // RK addition
#include <R_ext/BLAS.h> // RK addition
#include <math.h>*/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include "xxhash.h"
#define RATIO 10

// structs 
typedef struct {
    int x;
    int y;
} pair;

typedef struct {
	double * value_array; // log perms
	pair * array; // array with index
	int * table;
	size_t used_len; // current size of dict / number of vals
	size_t table_size; // current max size
	//SEXP arraySEXP;
	//SEXP value_arraySEXP;
} dictionary;


// Function declaration 

size_t hash(pair, int);
double Clog_sum_exp(double *, int, double);
void print_matrix(int, int , double * );
void print_int_vector(int ,  int * );
void print_float_vector(int ,  double * );

SEXP C_get_log_permanents(SEXP , SEXP , SEXP , SEXP , SEXP,SEXP);

int nonzero_perm(double * , double * , double * , int );

void get_union(int, double *, double * , int * , double *);


void update_dict(pair, double, dictionary *);
void print_sparse_matrix(dictionary * );

void sparse_reverse_tt(dictionary * ,  dictionary *,
						int, int,
						int * , int * , int *, double *, int, int *,
						int * );

void sparse_reverse_bt(dictionary * ,  dictionary * ,
						int , int ,
						int * , int *, int *, double *, int , int * ,
						int * );

void sparse_reverse_lm(dictionary *,  dictionary * ,
						int , int,
						int * , int * , int * , double * , int , int * ,
						int * );

void sparse_reverse_bs(dictionary * ,  dictionary * ,
						int, int ,
						int * , int *, int * , double * , int , int * ,
						int * );
void sparse_reverse_rm(dictionary * ,  dictionary * ,
						int, int ,
						int * , int *, int * , double * , int , int * ,
						int * );
void sparse_reverse_ts(dictionary * ,  dictionary * ,
						int, int ,
						int * , int *, int * , double * , int , int * ,
						int * );
void sparse_get_reduced_log_subperms(dictionary *, int * , int * , int *, 
				double * , int, int *, int *);

dictionary* sparse_reverse_reduction(dictionary * , dictionary * , int * , 
					   int * , int * , int *, int , int *, int * ,
			           int * , int * , double * );

double Csparse_log_sum_exp(dictionary *);


double lookup_logperm(dictionary *, pair );

void add_to_dictionary(dictionary * , pair, double);
dictionary * init_dictionary(size_t);

void free_dictionary(dictionary *);

void expand_dictionary(dictionary * );
void nullset_dictionary(dictionary * );

// reverse top trim
void reverse_tt(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				int * );

// Reverse bottom trim
void reverse_bt(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				int * );

// Reverse bottom split
void reverse_bs(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				double * , int * );
// Reverse top split
void reverse_ts(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				double * , int * );


// Reverse left merge
void reverse_lm(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				double * , int * );
// Reverse right merge
void reverse_rm(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				double * , int * );

void get_reduced_log_subperms(double * , int *, int *, int * , 
				double * , int, int *, int * );

void reduction(int * , int *, int * , int *, int, int *, int *,
			   int * , int * , int);

double * reverse_reduction(double * , double * , int * , 
					   int * , int * , int *, int , int *, int * ,
			           int *, int * , double * , 
			           double *);

void get_alphabetagamma(double * , int , double * , double * , double * , int , int *, 
    int *, int *, int * , int *, int );

// #define MADCONST 1.4826
// 
// void internal_matmult(double *,
//                       double *,
//                       double *,
//                       int, int, int, int);
// 
// void internal_matmultrightT(double *,
//                             double *,
//                             int , int );
// 
// void internal_matmultleftT(double * ,
//                            double * ,
//                            int , int );
// 
// double * internal_power_method(double *, int, double,
//                                int , double * , double * ,int  );
//   
// 
// void CUSUM(double *,double * , int , int , int );
// void singleCUSUM(double *, double *, int , int , int , int);
// 
// void rescale_variance(double *, double* , int, int, double *);
// 
// // Sorting
//     
// void insertSort (double * , int , int );
// 
// 
// void sort_k_largest(double * , int , int, int);
// void sort_k_largest_abs(double * , int , int, int);
// void sort_k_largest_int(int*, int , int, int);
// SEXP sort_k_largest_R(SEXP, SEXP, SEXP, SEXP);
// 
// 
// 
// 
// 
// //function to swap variable
// void swap(double*, double*);
// int partition (double *, int, int);
// int quickselect(double * , int, int, int);
// void rec_partial_quicksort(double *, int, int, int);
// 
// void partial_quicksort(double *, int , int);
// 
// 
// SEXP partial_quicksort_R(SEXP , SEXP, SEXP);


