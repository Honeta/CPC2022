#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix_def.h"
#include "mpi.h"
#include "host.h"

#define MAX_ITERATIONS 200

// ldu_data
LduMatrix ldu_matrix;
// csr_data
CsrMatrix csr_matrix;
double *x; //系数向量
double *b; //结果向量
double *source;

int comm_size = 0;
int my_rank = -1;
double spmv_time_cost[MATRIX_NUM];

void main_calc(int matrix_index) {
    int mesh_num = MESH_NUM[matrix_index];
    int row_num = ROW_NUM[matrix_index];
    int data_num = DATA_NUM[matrix_index];
    INFO("Mesh num             :  %d\n", mesh_num);
    INFO("Mesh >> decompose mesh!\n");
    decompose_mesh(mesh_num, row_num, data_num, ldu_matrix);

    INFO("Init >> init vector!\n");
    init_vector(row_num, x, b, source);

    INFO("Matrix >> matrix tranform!\n");
    ldu_to_csr(ldu_matrix, csr_matrix);

    INFO("Iteration >> Main Iterations start!\n");
    
    double xb = 1.0;
    for(int iter = 0; iter < MAX_ITERATIONS; iter++){
        if(iter % 100 == 0) {
            INFO("Iteration >> [%d/%d]\n", iter, MAX_ITERATIONS);
        }

        double spmv_time_start = MPI_Wtime();
        spmv(csr_matrix, x, b, my_rank);
        double spmv_time_end = MPI_Wtime();
        spmv_time_cost[matrix_index] += spmv_time_end - spmv_time_start;

        update_x(row_num, source, b, x, xb);
        if( (iter + 1 ) % 100 == 0) {
            INFO("Iteration >> [%d/%d] done!\n", iter, MAX_ITERATIONS);
        }
    }

    INFO("Write result!\n");
    write_result(row_num, data_num, x);
    
    INFO("Check >> check result!\n");
    if(check_result(row_num, data_num, x)) {
        INFO("Check >> check result correct!\n");
    } else {
        INFO("Check >> check result uncorrect!\n");
    }
    INFO("Time: spmv time cost: %.4lfs\n\n", spmv_time_cost[matrix_index]);

    free_matrix(ldu_matrix, csr_matrix);
    free_vector(x, b, source);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    INFO("Number of Process    :  %d\n\n", comm_size);

    double total_time_start =  MPI_Wtime();
    for(int i = 0; i < MATRIX_NUM; i++){
        INFO("<<================ Matrix index  :  %d ================>>\n", i);
        main_calc(i);
    }
    double total_time_end =  MPI_Wtime();

    double total_spmv_time_cost = .0;
    for(int i = 0; i < MATRIX_NUM; i++){
        INFO("Time: matrix %d spmv time cost: %.4lfs\n", i, spmv_time_cost[i]);
        total_spmv_time_cost += spmv_time_cost[i];
    }
    INFO("Time: total spmv time cost: %.4lfs\n", total_spmv_time_cost);
    INFO("Time: total time cost: %.4lfs\n", total_time_end - total_time_start);

    return 0;
}