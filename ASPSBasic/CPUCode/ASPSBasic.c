/**
 * 	Summary:
 *
 *  i-th column of b is in fmem
    a is in lmem
	small loop because of s
 *
 */
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <float.h>

#include <MaxSLiCInterface.h>
#include "Maxfiles.h"
#include "common.h"

// show help
int help_flag = 0;

// number of rows / columns. Matrix size is n*n
int n = 8;

/* tracing
	0 - prints only n, sum of result, realtime, cputime
	1 - prints input, output, final result
	2 - tests correctness of result and prints if test passed
*/
int trace = 2;

// elements of matrix are in -range  to +range interval
float range = 100.0;

void help(const char * cmd) {
    printf("Usage: %s [filename]\n", cmd);
    printf("\nOptions:\n");
    printf("  -h, --help\n\tPrint short help\n");
    printf("  -n, --size\n\tSize n of matrix\n");
    printf("  -r, --range\n\tRange of elements\n");
    printf("  -t, --trace\n\tTrace level: 0,1,2\n");

};

struct option options[] = {
	{"help",	required_argument, 0, 'h'},
	{"size",	required_argument, 0, 'n'},
	{"trace",	required_argument, 0, 't'},
	{"range",	required_argument, 0, 'r'},
	{0,0,0,0}
};

#define SHORTOPT "hn:t:r:"

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
			case '?':
				error(1, "Invalid1 option '%c'", optopt);
			default:
				abort();
		}
	}
	if (help_flag) {
		help(argv[0]);
		exit(0);
	}
}

void min_plus(int n, float *a, float *b, float *res){
	float *temp = malloc(n*n* sizeof(float));

	for (int i=0; i<n*n; i++) {
		temp[i] = FLT_MAX;
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				if (temp[i * n + j]> a[i * n + k] + b[k * n +j])
					temp[i * n + j] = a[i * n + k] + b[k * n +j];
			}
	  	}
	}

	for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				res[i*n+j] = temp[i*n+j];
		}
	}

}


void print_matrix(int n, float *mat){
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			if (mat[i*n+j]<FLT_MAX-10)
				printf("%f " , mat[i*n+j]);
			else
				printf(" inf ");
		}
		printf("\n");
	}
}


void apsp(int n, float *mat, float *res) {

	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			res[i*n+j] = mat[i*n+j];
		}
	}

	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (res[i * n + k] + res[k*n + j] < res[i * n +j]) {
					res[i * n + j] = res[i * n + k] + res[k * n + j];
				}
			}
	  	}
		printf("\n");
	}

}



int main(int argc, char * argv[])
{
	parse_args(argc, argv);
	const int size = n * n;
	int dataSizeBytes = size * sizeof(float);

	float *mat = malloc(dataSizeBytes);

	float *matTrans = malloc(dataSizeBytes);
	float *vector;

	float *output = malloc(dataSizeBytes);
	float *outputTrans = malloc(dataSizeBytes);
	float *expected = malloc(dataSizeBytes);

	for (int i=0; i<n*n; i++) {
		mat[i] = i+1;
	}

	generateMatrix(n, mat, range);

	/*for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			if (i==j){
				mat[i*n+j] = 0;
			}
			else {
				mat[i*n+j] = i+1;
			}
		}
	}

	float f = FLT_MAX;
	float matC[] =
		 {0, 5, f, 2, f, 1, 5, 3,
		  f, 0, 2, f, f, 1, 2, 3,
		  3, f, 0, f, 7, 1, 2, 3,
		  f, f, 4, 0, 1, 1, 3, 3,
		  1, 3, f, f, 0, 1, 2, f,
		  1, 3, f, f, 5, 0, 2, f,
		  1, 3, f, f, 1, 1, 0, f,
		  1, 3, f, f, 10, 1, 2, 0,
		 };

	 */

	timing_t timer;
	timer_start(&timer);

/*	for (int i=0; i<n; ++i) {
		vector = &matTrans[n*i];
		MatMatMultiply(n, matC, vector, &output[n*i]);
	}
	transpose(n, output, outputTran
	);
*/

	transpose(n, mat, matTrans);
	int j = 1;
	//while (j < n-1) {
		printf("mat C\n");
		printMatrix(n, mat);
		printf("mat trans\n");
		printMatrix(n, matTrans);

		for (int i=0; i<n; ++i) {
			vector = &matTrans[n*i];
			MatMatMultiply(n, mat, vector, &output[n*i]);
		}
		transpose(n, output, mat);
		matTrans = output;
		j *= 2;
	//}

	transpose(n, output, outputTrans);

	timer_stop(&timer);
	float sum = sumMat(size, outputTrans);
	printf("%d %f %ld %ld\n", n, sum, timer.realtime, timer.cputime);

	int status = 0;

	if (trace == 1) {

		printf("\nMatrix in \n");
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				printf("%f " , mat[i*n+j]);
			}
			printf("\n");
		}

		printf("\n\nResult\n");
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				printf("%f " , output[i*n+j]);
			}
			printf("\n");
		}
	}
	else if (trace == 2) {
		apsp(n, mat, expected);
		printMatrix(n, outputTrans);
		int status = check(size, outputTrans, mat);
		if (status) {
			printf("Test failed.\n");
			status = 1;
		}
		else
			printf("Test passed OK!\n");
	}

	/*
	free(mat);
	free(matTrans);
	free(res);
	free(output);
	free(outputTrans);
	free(expected);
 */
	return status;
}

