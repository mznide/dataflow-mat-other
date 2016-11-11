/**
 *
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
int st_vec = 60;


// number of rows / columns. Matrix size is n*n
int n = 4;
#define ALIGN_BURST		384
#define C MatMatMultiply_C

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

};

struct option options[] = {
	{"help",	required_argument, 0, 'h'},
	{"size",	required_argument, 0, 'n'},
	{"trace",	required_argument, 0, 't'},
	{"range",	required_argument, 0, 'r'},
	{"num_of_vec",	required_argument, 0, 'v'},
	{0,0,0,0}
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
				st_vec = atoi(optarg);
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

void mulMatVec(int n, float *mat, float *vec, float *res){
	for (int i = 0; i < n; i++) {
		float sum = 0;
		for (int j = 0; j < n; j++)
			sum += mat[i * n + j] * vec[j];
		res[i] = sum;
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


void transform_together(int n, float *inputMatrix, float *matrixTransformed, int vectorSize, int C1){
	int count = 0;
	for (int yy = 0; yy < n; yy += C1) {
		for ( int x = 0; x < n/vectorSize; ++x) {
			for ( int y = yy; y < yy + C1; ++y) {
				for (int v = 0; v < vectorSize; v++){
					matrixTransformed[count] = inputMatrix[y * n + x*vectorSize+v];
					count++;
					//printf("yy %d x %d y %d v %d \n", yy, x,y,v);
				}
			}
		}
	}
	//printf("n %d count %d \n", n, count);
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

int gcd(int a, int b) {
	if (b == 0) {
		return a;
	}
	else {
		return gcd(b, a % b);
	}
}

int main(int argc, char * argv[])
{
	parse_args(argc, argv);
	const int size = n * n;

	//int n1 = calc_align(n, VECTOR_SIZE*C <ALIGN_BURST/sizeof(float)? ALIGN_BURST/sizeof(float) : VECTOR_SIZE*C);

	int hcf = gcd(VECTOR_SIZE, 96);
	int lcm = (VECTOR_SIZE*96)/hcf;
	int n1 = calc_align(n, lcm);

	int size1 = n1*n1;

	int data_size_bytes = size * sizeof(float);
	int data_size_bytes1 = size1 * sizeof(float);

	float *mat_a = malloc(data_size_bytes);

	generate_matrix(n, mat_a, range);

	timing_t timer;
	timer_start(&timer);

	float *mat_a_aligned = malloc(data_size_bytes1);
	float *mat_a_trans = malloc(data_size_bytes1);

	alignMatrix(n, n1, mat_a, mat_a_aligned);

	//float vektorji_in[st_vec][n1];
	//float vektorji_out[st_vec][n1];
	//float vektorji_expected[st_vec][n];

	float *vektorji_in = malloc(sizeof(float)*st_vec*n1);
	float *vektorji_out = malloc(sizeof(float)*st_vec*n1);
	float *vektorji_expected = malloc(sizeof(float)*st_vec*n1);



	for (int i=0; i<st_vec; i++){
		generate_vector(n, &vektorji_in[i*n], range);
	}

	transform_together(n1, mat_a_aligned, mat_a_trans, VECTOR_SIZE, C);

	max_file_t * maxfile = MatMatMultiply_init();
	max_engine_t * engine = max_load(maxfile, "*");

	// write LP to LMem
	MatMatMultiply_writeLMem_actions_t writeact;
	writeact.param_address = 0;
	writeact.param_nbytes = data_size_bytes1;
	writeact.instream_cpu_to_lmem = mat_a_trans;
	MatMatMultiply_writeLMem_run(engine, &writeact);

	MatMatMultiply_actions_t actions;
	actions.param_matrixLength = n1*n1;
	actions.param_n = n1;

	for (int i=0; i<st_vec; i++) {
		actions.instream_vectorInput = &vektorji_in[i*n];
		actions.outstream_output = &vektorji_out[i*n];
		MatMatMultiply_run(engine, &actions);
	}

	max_unload(engine);

	timer_stop(&timer);

	printf("%d %d %ld %ld\n", n, st_vec, timer.realtime, timer.cputime);

	int status = 0;

	if (trace == 1) {
		printf("\nMatrix A\n");
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				printf("%f " , mat_a[i*n+j]);
			}
			printf("\n");
		}

		printf("\nVejtorji in \n");
		for (int i=0; i<st_vec; i++){
			printf("Vektor %d ", i);
			for (int j=0; j<n; j++){
				printf("%f " , vektorji_in[i*n+j]);
			}
			printf("\n");
		}

		printf("\nVejtorji out \n");
		for (int i=0; i<st_vec; i++){
			printf("Vektor %d ", i);
			for (int j=0; j<n; j++){
				printf("%f " , vektorji_out[i*n+j]);
			}
			printf("\n");
		}
	}
	else if (trace == 2) {
		int status = 0;
		for (int i=0; i<st_vec; i++) {
			mulMatVec(n, mat_a, &vektorji_in[i*n], &vektorji_expected[i*n]);
		}
		for (int i=0; i<st_vec; i++) {

			status = check(n, &vektorji_out[i*n], &vektorji_expected[i*n]) || status;
			printf("\n \n");
		}

		if (status) {
			printf("Test failed.\n");
		}
		else
			printf("Test passed OK!\n");
	}

	free(mat_a);
	free(mat_a_trans);
	free(mat_a_aligned);
	free(vektorji_in);
	free(vektorji_out);
	free(vektorji_expected);

	return status;
}
