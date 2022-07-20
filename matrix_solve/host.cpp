#include "host.h"
#include "matrix_def.h"
#include <athread.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG(...) fprintf(stderr, "debug: " __VA_ARGS__)

// clang format

int auto_tuner_impl(int N, int NNZ, int L, int *D) {
  float x = (float)NNZ / N;
  float sum = 0;
  for (int l = 1; l <= L; l++) {
    sum = sum + D[l] * (l - x) * (l - x);
  }
  float sigma = sqrt(sum / N);
  int K = L;
  int *DD;
  DD = (int *)malloc((L+1) * sizeof(int));
  DD[0] = 0;
  for (int i = 2; i <= L; i++) {
    for (int j = 1; j <= L; j++) {
      DD[j] = D[j];
    }
    int count = 0;
    for (int j = i + 1; j <= L; j++) {
      DD[i] = DD[i] + j / i * D[j];
      if (j % i == 0) {
        count = count + (j / i - 1) * D[j];
      } else {
        DD[j % i] = DD[j % i] + D[j];
      }
    }
    float xx = NNZ / (N + count);
    float summ = 0;
    for (int l = 1; l <= i; l++) {
      summ = summ + DD[l] * (l - xx) * (l - xx);
    }
    float sigma2 = sqrt(summ / (N + count));
    if (sigma2 <= sigma /* && enough storage available*/) {
      sigma = sigma2;
      K = i;
    }
  }
  return K;
}

int auto_tuner(int n, int nnz, int *cr_idx) {

  int *cnt = (int *)calloc(n, sizeof(int));
  for (int i = 0; i < nnz; i++) {
    assert(cr_idx[i] < n);
    cnt[cr_idx[i]]++;
  }

  int max_len = 0;
  for (int i = 0; i < n; i++) {
    if (cnt[i] > max_len) {
      max_len = cnt[i];
    }
  }

  int *D = (int *)calloc(max_len + 1, sizeof(int));

  for (int i = 0; i < n; i++) {
    assert(cnt[i] <= max_len);
    D[cnt[i]]++;
  }

  int ret = auto_tuner_impl(n, nnz, max_len, D);

clear:
  free(cnt);
  free(D);

  return ret;
}

void csr2csc(int n, int nnz, double *r_vals, int *r_cols, int *r_offs,
             double *c_vals, int *c_rows, int *c_offs) {

  for (int i = 0; i < nnz; i++) {
    if(r_cols[i] >= n) {
        //???????????????????????????????????????????????/
        // COL_IDX overflow?
        DEBUG("cols: %d, r_col: %d\n", n, r_cols[i]);
        assert(r_cols[i] < n);
    }
    c_offs[r_cols[i]]++;
  }

  // cumsum the nnz per column to get c_offs[]
  for (int col = 0, cumsum = 0; col < n; col++) {
    int temp = c_offs[col];
    c_offs[col] = cumsum;
    cumsum += temp;
  }

  for (int row = 0; row < n; row++) {
    int start = r_offs[row];
    int end = row == n - 1 ? nnz : r_offs[row + 1];
    for (int jj = start; jj < end; jj++) {
      int col = r_cols[jj];
      int dest = c_offs[col];

      c_rows[dest] = row;
      c_vals[dest] = r_vals[jj];

      c_offs[col]++;
    }
  }

  for (int col = 0, last = 0; col <= n; col++) {
    int temp = c_offs[col];
    c_offs[col] = last;
    last = temp;
  }
}

void init_cm_struct() {}

void init_ra_struct() {}

void cm_process() {}

void ra_process() {}

void CAspmv(const CsrMatrix &csr_matrix, double *x, double *b) {
  // DEBUG("%d %d\n", csr_matrix.rows, csr_matrix.row_off[csr_matrix.rows]);
  double *c_vals = (double *)malloc(csr_matrix.data_size * sizeof(double));
  int *c_rows = (int *)malloc(csr_matrix.data_size * sizeof(int));
  int *c_offs = (int *)calloc(csr_matrix.rows, sizeof(int));

//   csr2csc(csr_matrix.rows, csr_matrix.data_size, csr_matrix.data,
//           csr_matrix.cols, csr_matrix.row_off, c_vals, c_rows, c_offs);

  int cm_tuned_size =
      auto_tuner(csr_matrix.rows, csr_matrix.data_size, csr_matrix.cols);

  init_cm_struct();

  cm_process();

  init_ra_struct();

  ra_process();

clear:
  free(c_vals);
  free(c_rows);
  free(c_offs);
}

void naive_spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
  // deprecated
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

void spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    CAspmv(csr_matrix, x, b);
}

