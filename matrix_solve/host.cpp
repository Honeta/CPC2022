#include <athread.h>
#include "matrix_def.h"
#include "host.h"

const int nslaves = 64;

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

int GetRowIdx(int n, int x, int *off) {
    int l = 0;
    int r = n-1;
    if(x >= off[n-1]) return n-1;
    while(l < r) {
        int mid = (l+r+1) >> 1;
        if(off[mid] > x) {
            r = mid - 1;
        } else {
            l = mid;
        }
    }
    return l;
}


void spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    // naive_spmv(csr_matrix, x, b);
    int nnz = csr_matrix.data_size;
    int n = csr_matrix.rows;
    int *cols = csr_matrix.cols;
    double *value = csr_matrix.data;
    int *row_off = csr_matrix.row_off;
    
    int load = nnz / nslaves;

    int *start = (int*)malloc(sizeof(int)*nslaves);
    int *end = (int*)malloc(sizeof(int)*nslaves);
    int *L = (int*)malloc(sizeof(int)*nslaves);
    int *R = (int*)malloc(sizeof(int)*nslaves);
    double *midans = (double*)malloc(sizeof(double)*nslaves*2);

    for(int i = 0; i < nslaves; i++) {
        L[i] = i * load;
    }
    for(int i = 1; i < nslaves; i++) {
        R[i-1] = start[i];
    }
    R[nslaves-1] = nnz;

    for(int i = 0; i < nslaves; i++) {
        start[i] = GetRowIdx(n, L[i], row_off);
    }   
    for(int i = 1; i < nslaves; i++) {
        end[i-1] = start[i];
    }
    end[nslaves-1] = n;

    // start & end: edge case?

    SlaveArgs args = SlaveArgs{};
    args.x = x;
    args.b = b;
    args.row_ptr = row_off;
    args.L = L;
    args.R = R;
    args.start = start;
    args.end = end;
    args.cols = cols;

    // start slave func




    // reduce ans
    for(int i = 1; i < nslaves; i++) {
        int idx = start[i];
        b[idx] = midans[i<<1] + midans[i<<1|1];
    }

    free(start);
    free(end);
    free(L);
    free(R);
    free(midans);
}
