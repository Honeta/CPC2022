#ifndef _HOST_H_
#define _HOST_H_

#include "matrix_def.h"

void decompose_mesh(int mesh_num, int row_num, int data_num, LduMatrix &ldu_matrix);
void init_vector(int row_num, double* &x, double* &b, double* &source);

void ldu_to_csr(const LduMatrix &ldu_matrix, CsrMatrix &csr_matrix);    //数据转换
void spmv(const CsrMatrix &csr_matrix, double *x, double *b);           //计算
void update_x(int row_num, double *source, double *b, double *x, double &xb);

int check_result(int row_num, int data_num, double *x);
void write_result(int row_num, int data_num, double *x);

void free_matrix(LduMatrix &ldu_matrix, CsrMatrix &csr_matrix);
void free_vector(double* x, double* b, double* source);

struct AlbusArgs {
  int *col_idx;
  int *row_ptr;
  double *mtx_val;
  double *vec_val;
  double *mid_ans;
  int *start;
  int *end;
  int *block_size;
  int thread_nums;
  double *mtx_ans;
};

#endif