/**
 * 	Summary:
 *     CPU code for Matrix*Matrix multiplication.
 *     This solution uses FMem.
 *     This solution uses offset equal to size of the matrix so it's not suitable for big matrixes.
 */
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>

#include <MaxSLiCInterface.h>
#include "Maxfiles.h"
#include "common.h"
#include "float.h"

const int TILE_SIZE = MatMatMultiply_tileSize;


// show help
int help_flag = 0;
float f = FLT_MAX;

// number of rows / columns. Matrix size is n*n
int n = 4;

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

static int numOfTiles(int n) {
	return (n + TILE_SIZE - 1) / TILE_SIZE;
}

void generate_matrix1(int n, float *mat, int range)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			mat[i * n + j] = rand_int(range);
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

//floyd-warshall
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

/*
 * right order for dataflow input
 */
void reorder(int n, float *A, float *reordered) {
	int nTiles = numOfTiles(n);

	int pos = 0;
	for (int i = 0; i < nTiles; ++i) {
		for (int j = 0; j < nTiles; ++j) {
			for (int x = 0; x < TILE_SIZE; ++x) {
				int row = i*TILE_SIZE + x;
				for (int y = 0; y < TILE_SIZE; ++y) {
					int col = j*TILE_SIZE + y;
					reordered[pos] = A[row*n+col];
					pos += 1;
				}
			}
		}
	}
}

int calc_align(int n, int align) {
	return n / align * align + (n % align > 0) * align;
}

void align_matrix(int n, int n_aligned, float *mat, float *mat_aligned)
{
	for (int i = 0; i < n_aligned; i++) {
		for (int j = 0; j < n_aligned; j++){
			if (i <n && j<n) {
				mat_aligned[i * n_aligned + j] = mat[i * n + j];
			}
			else {
				mat_aligned[i * n_aligned + j] = f;
			}
		}
	}
}

void reverse_align_matrix(int n, int n_aligned, float *mat_aligned, float *mat_final)
{
	for (int i = 0; i < n_aligned; i++) {
		for (int j = 0; j < n_aligned; j++){
			if (i <n && j<n) {
				mat_final[i * n + j] = mat_aligned[i * n_aligned + j];
			}
		}
	}
}

void tile_multiply(int n1, int nTiles, int TILE_SIZE, float *mat_a, float *mat_b, float *mat_c_aligned) {

	int dataSizeBytes1 = n1 * n1* sizeof(float);
	float *a_reordered = malloc(dataSizeBytes1);
	float *b_reordered = malloc(dataSizeBytes1);

	reorder(n1, mat_a, a_reordered);
	reorder(n1, mat_b, b_reordered);

	float *output_tile = malloc(TILE_SIZE*TILE_SIZE*sizeof(float));
	float *sumA = malloc(TILE_SIZE*TILE_SIZE*sizeof(float));


	for (int i=0; i<nTiles; i++) {
		for (int j=0; j<nTiles; j++) {

			for (int z=0; z<TILE_SIZE*TILE_SIZE; z++)
				sumA[z]=FLT_MAX;

			for (int k=0; k<nTiles; k++){
				int a_tile = (i*nTiles+k)*TILE_SIZE*TILE_SIZE;
				int b_tile = (k*nTiles+j)*TILE_SIZE*TILE_SIZE;

				MatMatMultiply(TILE_SIZE*TILE_SIZE, &a_reordered[a_tile], &b_reordered[b_tile], output_tile);

				for (int z=0; z<TILE_SIZE*TILE_SIZE; z++) {
					sumA[z] = sumA[z] < output_tile[z]? sumA[z] : output_tile[z];
				}
			}

			int zacetek = i*n1*TILE_SIZE+j*TILE_SIZE;
			for (int i1=0; i1<TILE_SIZE; i1++) {
				for (int j1=0; j1<TILE_SIZE; j1++) {
					mat_c_aligned[zacetek+i1*n1+j1] = sumA[i1*TILE_SIZE+j1];
				}
			}
		  }
	}


	free(sumA);
	free(a_reordered);
	free(b_reordered);
	free(output_tile);
}

int main(int argc, char * argv[])
{
	parse_args(argc, argv);

	const int size = n*n;
	int n1 = calc_align(n, TILE_SIZE);
	int dataSizeBytes = size * sizeof(float);
	int dataSizeBytes1 = n1 * n1* sizeof(float);

	float *mat_a = malloc(dataSizeBytes);

	float *mat_final = malloc(dataSizeBytes1);

	generate_matrix1(n, mat_a, 15);

	//nicle po diagonali
	for (int i=0; i<n; i++){
		mat_a[i*n+i] = 0;
	}

	timing_t timer;
	timer_start(&timer);

	//adds zeros to make matrix divisable with tile size
	float *mat_a_aligned = malloc(dataSizeBytes1);
	float *mat_c_aligned = malloc(dataSizeBytes1);

	float *expected = malloc(dataSizeBytes);

	align_matrix(n, n1, mat_a, mat_a_aligned);

	int nTiles = numOfTiles(n1);

	int i = 1;
	//i=1 do log2n
	while (i<=n) {
		tile_multiply(n1, nTiles, TILE_SIZE, mat_a_aligned, mat_a_aligned, mat_a_aligned);
		//tile_multiply(n1, nTiles, TILE_SIZE, mat_a_aligned, zacetna, mat_a_aligned);
		i *= 2;
		//i += 1;
		//printf("korak %d \n", i);
		//print_matrix(n, mat_final);
	}

	reverse_align_matrix(n, n1, mat_a_aligned, mat_final);
	timer_stop(&timer);

	float sum = sum_mat(size, mat_final);

	printf("%d %f %ld %ld\n", n, sum, timer.realtime, timer.cputime);

	int status = 0;

	if (trace == 1) {
		print_matrix(n, mat_a);
		print_matrix(n, mat_final);
	}
	else if (trace == 2) {
		apsp(n, mat_a, expected);
		int status = check(n*n, mat_final, expected);
		if (status) {
			printf("Test failed.\n");
			status = 1;
		}
		else
			printf("Test passed OK!\n");
	}

	free(mat_a);
	free(mat_a_aligned);
	free(mat_c_aligned);
	free(mat_final);
	free(expected);

	return status;
}
