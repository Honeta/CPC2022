#include "host.h"
#include "matrix_def.h"
#include <athread.h>

const int thread_num = 64;
extern "C" {
    void slave_func_albus_part1(void *);
    void slave_func_albus_part2(void *);
}

void naive_spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    // 实现样例
    for (int i = 0; i < csr_matrix.rows; i++) {
      int start = csr_matrix.row_off[i];
      int num = csr_matrix.row_off[i + 1] - csr_matrix.row_off[i];
      double result = .0;
      for (int j = 0; j < num; j++) {
        result += x[csr_matrix.cols[start + j]] * csr_matrix.data[start + j];
      }
      b[i] = result;
    }
}

int GetRowIdx(int n, int x, int *off) {
    int l = 0;
    int r = n - 1;
    if (x >= off[n - 1])
      return n - 1;
    while (l < r) {
      int mid = (l + r + 1) >> 1;
      if (off[mid] > x) {
        r = mid - 1;
      } else {
        l = mid;
      }
    }
    return l;
}

void spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    // naive_spmv(csr_matrix, x, b);
    int row_num = csr_matrix.rows;
    int nnz = csr_matrix.row_off[row_num];
    int *col_idx = csr_matrix.cols;
    double *mtx_val = csr_matrix.data;
    int *row_ptr = csr_matrix.row_off;

    int *start = (int *)malloc(sizeof(int) * thread_num);
    int *end = (int *)malloc(sizeof(int) * (thread_num + 1));
    int *block_size = (int *)malloc(sizeof(int) * (thread_num + 1));
    double *mid_ans = (double *)malloc(sizeof(double) * thread_num * 2);

    // Algorithm 4 : ALBUS load balancing settings
    block_size[thread_num] = nnz;
    start[0] = 0;
    end[thread_num - 1] = row_num;
    int t = nnz / thread_num;
    for (int i = 0; i < thread_num; i++) {
      block_size[i] = i * t;
      start[i] = GetRowIdx(row_num, block_size[i], row_ptr);
      end[i - 1] = start[i];
    }

    // start & end: edge case?

    // start slave func
    AlbusArgs args = AlbusArgs{};
    args.col_idx = col_idx;
    args.row_ptr = row_ptr;
    args.mtx_val = mtx_val;
    args.vec_val = x;
    args.mid_ans = mid_ans;
    args.start = start;
    args.end = end;
    args.block_size = block_size;
    args.thread_nums = thread_num;
    args.mtx_ans = b;

    CRTS_init();
    CRTS_athread_spawn(slave_func_albus_part1, &args);
    CRTS_athread_join();
    CRTS_athread_spawn(slave_func_albus_part2, &args);
    CRTS_athread_join();
    CRTS_athread_halt();
    // reduce ans
    // for(int i = 1; i < thread_num; i++) {
    //     int idx = start[i];
    //     b[idx] = midans[i<<1] + midans[i<<1|1];
    // }

    free(start);
    free(end);
    free(block_size);
    free(mid_ans);
}
