#include <athread.h>
#include "matrix_def.h"
#include "host.h"

struct Para {
    double *vec;
    double *mat_val;
    double *res;
    int *cols;
    int *block;
    int *r_ptr;
};

extern "C" {
    void slave_func_calc(Para *x);
}

void naive_spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    // 实现样例
    for(int i = 0; i < csr_matrix.rows; i++) {
        int start = csr_matrix.row_off[i];
        int num = csr_matrix.row_off[i+1] - csr_matrix.row_off[i];
        double result = .0;
        for(int j = 0; j < num; j++) {
            result += x[csr_matrix.cols[start+j]] * csr_matrix.data[start+j];
        }
        b[i]=result;
    }
}


void spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    int n = csr_matrix.rows;
    int nnz = csr_matrix.row_off[n];
    double *val = csr_matrix.data;
    int *cols = csr_matrix.cols;
    int *r_ptr = csr_matrix.row_off;

    const int nt = 64;

    int load = n / nt;
    int *block = (int*)malloc(sizeof(int)*(nt+1));
    for(int i = 0; i < nt; i++) {
        block[i] = load * i;
    }
    block[nt] = n;

    Para arg = Para{};
    arg.vec = x;
    arg.mat_val = val;
    arg.res = b;
    arg.cols = cols;
    arg.block = block;
    arg.r_ptr = r_ptr;
    CRTS_init();
    CRTS_athread_spawn(slave_func_calc, &arg);
    CRTS_athread_join();
    CRTS_athread_halt();

    free(block);

}
