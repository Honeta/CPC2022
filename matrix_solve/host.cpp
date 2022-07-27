#include <athread.h>
#include "matrix_def.h"
#include "host.h"


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
    // naive_spmv(csr_matrix, x, b);
    
}
