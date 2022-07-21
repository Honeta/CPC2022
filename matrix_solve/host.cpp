#include <athread.h>
#include <vector>
#include "matrix_def.h"
#include "host.h"

void block_sub_mat(const CsrMatrix &csr_matrix) {
  const int div = 8;
  const int tile = csr_matrix.rows / div;
  auto cnt = std::vector<std::vector<int>>(div, std::vector<int>(div, 0));
  for(int i = 0; i < csr_matrix.rows; i++) {
    int r = i / tile;
    int start = csr_matrix.row_off[i];
    int end = i+1 == csr_matrix.rows ? csr_matrix.data_size : csr_matrix.row_off[i+1];
    for(int j = start; j < end; j++) {
      if(csr_matrix.cols[j] > csr_matrix.rows) continue;
      int c = csr_matrix.cols[j] / tile;
      cnt[r][c]++;
    }
  }
  int ave = csr_matrix.data_size / div / div;
  for(int i = 0; i < div; i++) {
    for(int j = 0; j < div; j++) {
      printf("%d ", cnt[i][j] - ave);
    }
    printf("\n");
  }

  printf("\n");
}

void block_rows(const CsrMatrix &csr_matrix) {
  int nnz = csr_matrix.data_size;
  int n = csr_matrix.rows;
  // cols == rows
  auto cnt = std::vector<int>(64, 0);
  const int div = 64;
  int tile = n / div;
  int ave = nnz / div;
  for(int i = 0; i < div; i++) {
    int start = i * tile;
    int end = start + tile;
    if(end > n) {
      end = n;
    }
    for(int j = start; j < end; j++) {
      int ele_end = j + 1 == n ? nnz : csr_matrix.row_off[j+1];
      cnt[i] += ele_end - csr_matrix.row_off[j];
    }
  }
  for(int i = 0; i < div; i++) {
    printf("%d ", cnt[i]);
  }
  printf("\n");
}


void spmv(const CsrMatrix &csr_matrix, double *x, double *b, int rank) {
    // 实现样例
    if(rank != 0) return;
    if(csr_matrix.row_off[csr_matrix.rows] != csr_matrix.data_size) {
      INFO("overflow row_off");
    }
    
    // block_sub_mat(csr_matrix);

    // block_rows(csr_matrix);

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
