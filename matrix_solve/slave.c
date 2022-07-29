#include <slave.h>
#include "matrix_def.h"

struct Para {
    double *vec;
    double *mat_val;
    double *res;
    int *cols;
    int *block;
    int *r_ptr;
};


void slave_func_calc(struct Para *x) {
    int tid = CRTS_tid;
    int start = x->block[tid];
    int end = x->block[tid+1];
    for(int i = start; i < end; i++) {
        int l = x->r_ptr[i];
        int r = x->r_ptr[i+1];
        double res = .0;
        for(int j = l; j < r; j++) {
            res += x->vec[x->cols[j]] * x->mat_val[j];
        }
        x->res[i] = res;
    }
}