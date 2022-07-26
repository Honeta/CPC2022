#include <slave.h>
#include <crts.h>
#include "matrix_def.h"

#define ldm_size 256*1024

struct Csc_matrix {
    int n;
    int *col_ptr;       //index of the first element in a column
    int *col_emt;       //number of element in every column
    int *rows;          //row number of every non-0 element
    double *value;      //value of every non-0 element
    int data_size;

    int inc_cnt;
    int *splited_col_ptr;
};

typedef struct {
    struct Csc_matrix csc_matrix;
    int *P;
    double *splited_x;
    int Inc;
} CmPara;

struct Result_matrix {
    int n;
    int *row_ptr;
    int *row_emt;
    int *columns;
    double *value;
    int data_size;

    int inc_cnt;
    int *splited_row_ptr;
};

typedef struct {
    struct Result_matrix result_matrix;
    int *P;
    double *splited_b;
    int Inc;
} RaPara;

__thread_local crts_rply_t dma_rply = 0;
__thread_local unsigned int D_COUNT = 0;
// __thread_local double

// checkpoint
void CP(const void *a) { printf("%s\n", a); }

void slave_fun_cm(void *para) {
    CmPara cmpara;
    //split rules
    int *P;
    //start column and end column number
    int start, end;

    //get CPE id
    int tid = CRTS_tid;

    CP("E");

    //get para structure
    CRTS_dma_iget(&cmpara, para, sizeof(CmPara), &dma_rply);
    ++D_COUNT;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    //set task of every CPE
    P = cmpara.P;
    start = P[tid];
    end = P[tid + 1];
    
    //set ldm rounds
    int Inc = cmpara.Inc;
    //ldm rounds number of CPE
    int rounds = (end - start) / Inc;
    //remainder rows start number after rounds ldm round has processed
    int left = Inc * rounds + start;

    CP("E 0");
    //init x and col_ptr(splited) in ldm
    double * x = (double *)ldm_malloc(Inc * sizeof(double));
    int * col_ptr  = (int *)ldm_malloc((Inc + 1) * sizeof(int));
    double * data;
    int len;
    
    CP("E 1");
    for (int round = 0; round < rounds; ++round) {
        int round_col = start + round * Inc;
        //load x and col_ptr of this round to ldm
        CRTS_dma_iget(x, cmpara.splited_x + round_col, 
            Inc * sizeof(double), &dma_rply);
        CRTS_dma_iget(col_ptr, cmpara.csc_matrix.splited_col_ptr + round_col, 
            (Inc + 1) * sizeof(int), &dma_rply);
        D_COUNT += 2;
        CRTS_dma_wait_value(&dma_rply, D_COUNT);
        
        CP("E 2");
    
        //init data in ldm(data dependence:col_ptr)
        //! if error, test to pin length
        int start_idx = col_ptr[0];
        len = col_ptr[Inc] - start_idx;
        data = (double *)ldm_malloc(len * sizeof(double));

        //load value of this round to ldm
        CRTS_dma_iget(data, cmpara.csc_matrix.value + start_idx,
                len * sizeof(double), &dma_rply);
        ++D_COUNT;
        CRTS_dma_wait_value(&dma_rply, D_COUNT);

        CP("E 3");

        //multiply processing
        int idx = 0;
        for (int i = 0; i < Inc; ++i) {
            int num = col_ptr[i + 1] - col_ptr[i];
            for (int j = 0; j < num; ++j) {
                data[idx] = data[idx] * x[i];
                ++idx;
            }
        }
        //! try to using double buffer here(r&w on the same time)
        CRTS_dma_iput(cmpara.csc_matrix.value + start_idx, data, 
                len * sizeof(double), &dma_rply);
        ++D_COUNT;
        CRTS_dma_wait_value(&dma_rply, D_COUNT);
        ldm_free(data, len * sizeof(double));
    }

    //remainder processing
    CRTS_dma_iget(x, cmpara.splited_x + left, 
        (end - left) * sizeof(double), &dma_rply);
    CRTS_dma_iget(col_ptr, cmpara.csc_matrix.splited_col_ptr + left, 
        (end - left + 1) * sizeof(int), &dma_rply);
    D_COUNT += 2;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);
    len = col_ptr[end - left] - col_ptr[0];
    data = (double *)ldm_malloc(len * sizeof(double));
    CRTS_dma_iget(data, cmpara.csc_matrix.value + col_ptr[0],
            len * sizeof(double), &dma_rply);
    ++D_COUNT;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    CP("E 4");
    int idx = 0;
    for (int i = 0; i < end - left; ++i) {
        int num = col_ptr[i + 1] - col_ptr[i];
        for (int j = 0; j < num; ++j) {
            data[idx] = data[idx] * x[i];
            ++idx;
        }
    }
    CP("E 5");
    CRTS_dma_iput(cmpara.csc_matrix.value + col_ptr[0],
            data, len * sizeof(double), &dma_rply);
    ++D_COUNT;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);
    ldm_free(data, len * sizeof(double));
    ldm_free(x, Inc * sizeof(double));
    ldm_free(col_ptr, (Inc + 1) * sizeof(int));
    CP("E 6");
}

void slave_fun_ra(void *para) {
  RaPara rapara;
  int *P;
  int start, end;
  int tid = CRTS_tid;

  // get para structure
  CRTS_dma_iget(&rapara, para, sizeof(RaPara), &dma_rply);
  ++D_COUNT;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  // set task
  P = rapara.P;
  start = P[tid];
  end = P[tid + 1];

  int Inc = rapara.Inc;
  int rounds = (end - start) / Inc;
  int left = Inc * rounds + start;

  // init b and row_ptr(splited) in ldm
  double *b = (double *)ldm_malloc(Inc * sizeof(double));
  int *row_ptr = (int *)ldm_malloc((Inc + 1) * sizeof(int));
  double *data;
  int len, round_row;
  for (int round = 0; round < rounds; ++round) {
    // start row of this rounds
    round_row = start + round * Inc;
    CRTS_dma_iget(row_ptr, rapara.result_matrix.splited_row_ptr + round_row,
                  (Inc + 1) * sizeof(int), &dma_rply);
    ++D_COUNT;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    // init data in ldm(data dependence:row_ptr)
    //! if error, test to pin length
    int start_idx = row_ptr[0];
    len = row_ptr[Inc] - start_idx;
    data = (double *)ldm_malloc(len * sizeof(double));

    // load value of this round to ldm
    CRTS_dma_iget(data, rapara.result_matrix.value + start_idx,
                  len * sizeof(double), &dma_rply);
    ++D_COUNT;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    // add processing
    int idx = 0;
    for (int i = 0; i < Inc; ++i) {
      int num = row_ptr[i + 1] - row_ptr[i];
      b[i] = 0;
      for (int j = 0; j < num; ++j) {
        b[i] += data[idx];
        ++idx;
      }
    }
    //! try to using double buffer here(r&w on the same time)
    CRTS_dma_iput(rapara.splited_b + round_row, b, Inc * sizeof(double),
                  &dma_rply);
    ++D_COUNT;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);
    ldm_free(data, len * sizeof(double));
  }

  // remainder processing
  CRTS_dma_iget(row_ptr, rapara.result_matrix.splited_row_ptr + left,
                (end - left + 1) * sizeof(int), &dma_rply);
  ++D_COUNT;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);
  len = row_ptr[end - left] - row_ptr[0];
  data = (double *)ldm_malloc(len * sizeof(double));
  CRTS_dma_iget(data, rapara.result_matrix.value + row_ptr[0],
                len * sizeof(double), &dma_rply);
  ++D_COUNT;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  int idx = 0;
  for (int i = 0; i < end - left; ++i) {
    int num = row_ptr[i + 1] - row_ptr[i];
    b[i] = 0;
    for (int j = 0; j < num; ++j) {
      b[i] += data[idx];
      ++idx;
    }
  }
  CRTS_dma_iput(rapara.splited_b + left, b, Inc * sizeof(double), &dma_rply);
  ++D_COUNT;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);
  ldm_free(data, len * sizeof(double));
  ldm_free(b, Inc * sizeof(double));
  ldm_free(row_ptr, (Inc + 1) * sizeof(int));
}
