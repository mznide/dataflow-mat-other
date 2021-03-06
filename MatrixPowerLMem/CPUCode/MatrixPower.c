/**
 * V bistvu tam bo sam pršparala eno iteracijo branja in pisanja v LMem kar je bolj bogo...
 *
 * 	Summary:
 */
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>

#include <MaxSLiCInterface.h>
#include "Maxfiles.h"
#include "common.h"

#define VECTOR_SIZE MatMatMultiply_vectorSize

//TODO:  matrixes < 32

// show help
int help_flag = 0;

// number of rows / columns. Matrix size is n*n
int n = 4;

//power;
int power = 2;


#define ALIGN_BURST		384
/* tracing
	0 - prints only n, sum of result, realtime, cputime
	1 - prints input, output, final result
	2 - tests correctness of result and prints if test passed
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
    printf("  -t, --trace\n\tTrace level: 0,1,2\n");
    printf("  -p, --power\n\tMatrix power: 0,1,...\n");

};

struct option options[] = {
	{"help",	required_argument, 0, 'h'},
	{"size",	required_argument, 0, 'n'},
	{"trace",	required_argument, 0, 't'},
	{"range",	required_argument, 0, 'r'},
	{"power",	required_argument, 0, 'r'},
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
			case 'p':
				power = atoi(optarg);
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


//transforms matrix for vectorized multiplication
void transform (int n, float *input, float *matrixTransformed, int vectorSize){
	for (int v = 0; v < n; v=v+vectorSize) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < vectorSize; j++) {
				matrixTransformed[i*vectorSize +j + v*n] = input[j+i*n+v];

			}
		}
	}
}

int calc_align(int n, int align) {
	return n / align * align + (n % align > 0) * align;
}

void alignMatrix(int n, int n_aligned, float *mat, float *mat_aligned)
{
	for (int i = 0; i < n_aligned; i++) {
		for (int j = 0; j < n_aligned; j++){
			if (i <n && j<n) {
				mat_aligned[i * n_aligned + j] = mat[i * n + j];
			}
			else {
				mat_aligned[i * n_aligned + j] = 0;
			}
		}
	}
}


//res variable independent - can be same as input

void mulMatMat(int n, float *mat, float *matB, float *res){
	float *temp = malloc(n*n* sizeof(float));

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			float sum = 0;
			for (int k = 0; k < n; k++) {
				sum +=  mat[i*n+k] * matB[k*n+j];
			}
			temp[i*n+j] = sum;
	  }
	}
	for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				res[i*n+j] = temp[i*n+j];
		}
	}
}

void mat_power(int n, float *mat, int pow, float *res) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++){
			res[i*n+j] = i==j;
		}
	}

	while (pow>0) {
        if (pow % 2 == 0) {
            mulMatMat(n, mat, mat, mat);
            pow /= 2;
        }
        else {
        	mulMatMat(n, res, mat, res);
        	pow -= 1;
        }
	}

}


void reverseAlignMatrix(int n, int n_aligned, float *mat_aligned, float *mat_final)
{
	for (int i = 0; i < n_aligned; i++) {
		for (int j = 0; j < n_aligned; j++){
			if (i <n && j<n) {
				mat_final[i * n + j] = mat_aligned[i * n_aligned + j];
			}
		}
	}
}



int main(int argc, char * argv[])
{
	parse_args(argc, argv);
	const int size = n * n;
	int n1 = calc_align(n, ALIGN_BURST/sizeof(float));
	int size1 = n1*n1;

	int sizeBytes = size * sizeof(float);
	int sizeBytes1 = size1 * sizeof(float);

	float *mat = malloc(sizeBytes);

	float *matAligned = malloc(sizeBytes1);

	float *matTransformed = malloc(sizeBytes1);
	float *matTransposed = malloc(sizeBytes1);
	float *mat_final = malloc(sizeBytes);


	float *vector;
	float *output = malloc(sizeBytes1);

	float *expected = malloc(sizeBytes);

	generateMatrix(n, mat, range);
	for (int i=0; i<n*n; i++) {
			mat[i] = i+1;
	}
	alignMatrix(n, n1, mat, matAligned);


	timing_t timer;
	timer_start(&timer);


	int first = 1;
	float *res = malloc(sizeBytes1);

	int pow1 = power;
	while (pow1>0) {
		if (pow1 % 2 == 0) {
			transform(n1, matAligned, matTransformed, VECTOR_SIZE);
			transpose(n1, matAligned, matTransposed);

			MatMatMultiply_writeLMem(0, sizeBytes1, matTransformed);

			for (int i=0; i<n1; i++) {
				vector = &matTransposed[n1*i];
				MatMatMultiply(i, size1, n1, vector);
			}

			MatMatMultiply_readLMem(sizeBytes1, sizeBytes1, output);

			transpose(n1, output, matAligned);

			pow1 /=2;

		}
		else {
			printf("pow2\n");
			if (first) {
				for (int i=0; i<n1; i++) {
					for (int j=0; j<n1; j++) {
						res[i*n1+j] = matAligned[i*n1+j];
					}
				}
				first = 0;
			}
			else {
				//recimo tuki bo vedno res druga
				transform(n1, matAligned, matTransformed, VECTOR_SIZE);
				transpose(n1, res, matTransposed);

				MatMatMultiply_writeLMem(0, sizeBytes1, matTransformed);

				for (int i=0; i<n1; i++) {
					vector = &matTransposed[n1*i];
					MatMatMultiply(i, size1, n1, vector);
				}
				MatMatMultiply_readLMem(sizeBytes1, sizeBytes1, output);
				transpose(n1, output, res);
			}
			pow1 -= 1;
		}
	}

	reverseAlignMatrix(n, n1, res, mat_final);

	timer_stop(&timer);

	float sum = sumMat(size, output);
	printf("%d %f %ld %ld\n", n, sum, timer.realtime, timer.cputime);

	int status = 0;

	if (trace == 1) {
		printf("\nMatrix A\n");
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
		mat_power(n, mat, power, expected);
		printf("power %d", power);
		printMatrix(n, mat);
		printMatrix(n, expected);
		int status = check(size, mat_final, expected);
		if (status) {
			printf("Test failed.\n");
			status = 1;
		}
		else
			printf("Test passed OK!\n");

	}

	free(mat);
	free(matTransformed);
	free(matTransposed);
	free(matAligned);
	free(output);
	free(mat_final);
	free(expected);

	return status;
}
