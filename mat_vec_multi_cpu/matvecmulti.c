/**
 * 	Summary:
 *   CPU code for basic Matrix*vector multiplication.
 *	 matrix is read as a stream
 * 	 vector is stored in fmem
 *	 sum creates a small loop to wait for addition
 *   number of ticks is n*n*13 (13 because of adder component)
 */
 #define PLAINTIMER

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include "common.h"


// **************** arguments

void help(const char * cmd) {
    printf("Usage: %s [filename]\n", cmd);
    printf("\nOptions:\n");
    printf("  -h, --help\n\tPrint short help\n");
    printf("  -n, --size\n\tSize n of matrix\n");
    printf("  -r, --range\n\tRange of elements\n");
    printf("  -t, --trace \n\tTrace level: 0,1,2\n");

};


// show help
int help_flag = 0;

// number of rows / columns. Matrix size is n*n
int n = 4;

/* tracing
	0 - prints only n, sum of result, realtime, cputime
	1 - prints input, output, final result
*/
int trace = 0;

// elements of matrix are in -range  to +range interval
float range = 100.0;
int num_of_vec = 10;

struct option options[] = {
	{"help",	required_argument, 0, 'h'},
	{"size",	required_argument, 0, 'n'},
	{"trace",	required_argument, 0, 't'},
	{"range",	required_argument, 0, 'r'},
	{"num_of_vec",	required_argument, 0, 'v'},
};

#define SHORTOPT "hn:t:r:v:"

void parse_args(int argc, char * argv[]) {
	while (1) {
		int option_index = 0;
		int opt = getopt_long(argc, argv, SHORTOPT, options, &option_index);
		if (opt == -1) break;

		switch (opt) {
			case 'h':
				help_flag = 1;
				break;
			case 'n':
				n = atoi(optarg);
				break;
			case 't':
				trace = atoi(optarg);
				break;
			case 'r':
				range = atoi(optarg);
				break;
			case 'v':
				num_of_vec = atoi(optarg);
				break;
			case '?':
				error(1, "Invalid option '%c'", optopt);
			default:
				abort();
		}
	}
	if (help_flag) {
		help(argv[0]);
		exit(0);
	}
}


// **************** matrices


void generateMat(int n, float *mat, int range) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			mat[i * n + j] = rand_signdouble(range);
}


void generateVec(int n, float *vec, int range) {
	for (int i = 0; i < n; i++) 
		vec[i] = rand_signdouble(range);	
		
}


void printVec(int n, float *vec) {
	for (int i = 0; i < n; i++)
		printf("%.2f ", vec[i]);
	puts("\n");
}


void mulMatVec(int n, float *mat, float *vec, float *res){
	for (int i = 0; i < n; i++) {
		float sum = 0;
		for (int j = 0; j < n; j++)
			sum += mat[i * n + j] * vec[j];
		res[i] = sum;
	}
}


float sumMatDiag(int n, float *mat) {
	// trace of a matrix
	float sum = 0;
	for (int i = 0; i < n; i++)
		sum += mat[i * n + i];
	return sum;
}


float sumVec(int n, float *vec) {
	float sum = 0;
	for (int i = 0; i < n; i++)
		sum += vec[i];
	return sum;
}


int check(float *output, float *expected) {
	int status = 0;
	for (int i = 0; i < n; i++) {
		if (output[i] != expected[i]) {
			fprintf(stderr, "[%d] error, output: %f != expected: %f\n",
				i, output[i], expected[i]);
			status = 1;
		}
	}
	return status;
}


int main(int argc, char * argv[]) {
	parse_args(argc, argv);

	int matSizeBytes = n * n * sizeof(float);
	int vecSizeBytes = n * sizeof(float);

	float *mat = malloc(matSizeBytes);
	float *res = malloc(vecSizeBytes);

	//float (*sample) [num_of_vec];
	//float *s = malloc(sizeof(*sample) * n*2); 

	//float vektorji_in[num_of_vec][n];
	//float vektorji_out[num_of_vec][n];

	float *vektorji_in = malloc(n*num_of_vec*sizeof(float));
	float *vektorji_out = malloc(n*num_of_vec*sizeof(float));


	generateMat(n, mat, range);


	for(int i=0; i<num_of_vec; i++){
		generateVec(n, &vektorji_in[i*n], range);		
	}
	

	timing_t timer;
	timer_start(&timer);
	
	for(int i=0; i<num_of_vec;i++) {
		mulMatVec(n, mat, &vektorji_in[i*n], &vektorji_out[i*n]);
	}


	timer_stop(&timer);

	printf("%d %d %ld %ld\n", n, num_of_vec, timer.realtime, timer.cputime);

	if (trace == 1) {
		printf("\nMatrix\n");
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				printf("%f " , mat[i*n+j]);
			}
			printf("\n");
		}
		printf("\nVector \n");
		for (int i=0; i<n; i++){
			printVec(n, &vektorji_in[i*n]);
		}
		printf("\n\nResult\n");
		for (int i=0; i<num_of_vec; i++){
			printVec(n, &vektorji_out[i*n]);
		}
		printf("\n");
	}	
/*
	free(mat);
	free(vec);
	free(res);
*/
	return 0;
}
