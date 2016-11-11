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
#include <float.h>

// show help
int help_flag = 0;

// number of rows / columns. Matrix size is n*n
int n = 8;

float f = FLT_MAX;


/* tracing
	0 - prints only n, sum of result, realtime, cputime
	1 - prints input, output, final result
*/
int trace = 0;

// elements of matrix are in -range  to +range interval
float range = 100.0;


void help(const char * cmd) {
    printf("Usage: %s [filename]\n", cmd);
    printf("\nOptions:\n");
    printf("  -h, --help\n\tPrint short help\n");
    printf("  -n, --size\n\tSize n of matrix\n");
    printf("  -r, --range\n\tRange of elements\n");
    printf("  -t, --trace\n\tTrace level: 0,1\n");
    printf("  -p, --power\n\tMatrix power: 0,1,...\n");
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


// **************** matrices



void generateMat(int n, float *mat, int range) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			mat[i * n + j] = rand_int(range);
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

void minPlus(int n, float *mat, float *res){
	float *temp = malloc(n*n* sizeof(float));


	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			float min = f;
			for (int k = 0; k < n; k++) {
				if (mat[i*n+k]+mat[k*n+j]<min)
					min = mat[i*n+k]+mat[k*n+j];
			
			}
			temp[i*n+j] = min;
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
				printf("  %.0f " , mat[i*n+j]);
			else
				printf(" inf ");
		}
		printf("\n");
	}
	printf("\n");
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

float sumMat(int size, float *mat) {
	float sum = 0;
	for (int i = 0; i < size; i++)
		sum += mat[i];
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
	}

}



int main(int argc, char * argv[]) {
	parse_args(argc, argv);
	int size = n*n;
	int matSizeBytes = size * sizeof(float);

	float *mat = malloc(matSizeBytes);
	float *res = malloc(matSizeBytes);
	float mat1[] =
	 {0, 5, f, 2, f, 1, 5, 3,
	  f, 0, 2, f, f, 1, 2, 3,
	  3, f, 0, f, 7, 1, 2, 3,
	  f, f, 4, 0, 1, 1, 3, 3,
	  1, 3, f, f, 0, 1, 2, f,
	  1, 3, f, f, 5, 0, 2, f,
	  1, 3, f, f, 1, 1, 0, f,
	  1, 3, f, f, 10, 1, 2, 0,
	 };
	
	generateMat(n, mat, 8);


	if (n != 8)
		generateMat(n, mat, range);
	//nicle
	for(int i=0;i<n;i++){
		mat[i*n+i]=0;	
	}
	timing_t timer;
	timer_start(&timer);
	//multMatMat(n, matA, matB, res);

	int pow = n;
	for (int i=1; i<=pow; i=i*2) {
		minPlus(n, mat, mat);
	}

	timer_stop(&timer);

	float tr = sumMat(size, res);
	printf("%d %.2f %ld %ld\n", n, tr, timer.realtime, timer.cputime);
	
	if (trace == 1) {
		printf("\nMatrix in\n");
		print_matrix(n, mat);

		printf("\n\nResult\n");
		print_matrix(n, mat1);
	}
	free(res);
	free(mat);
	return 0;
}
