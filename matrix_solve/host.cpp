#include <athread.h>
#include "matrix_def.h"
#include "host.h"

// x : vector
// b : result
void spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    // 实现样例
    for(int i = 0; i < csr_matrix.rows; i++) {
        //获取该行第一个非0元 在所有非0元中的索引
        int start = csr_matrix.row_off[i];
        //获取该行非0元数量
        int num = csr_matrix.row_off[i+1] - csr_matrix.row_off[i];
        double result = .0;
        for(int j = 0; j < num; j++) {
            result += x[csr_matrix.cols[start+j]] * csr_matrix.data[start+j];
        }
        b[i]=result;
    }
}
