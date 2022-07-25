#include <slave.h>
#include <crts.h>
#include "matrix_def.h"

#define ldm_size 256*1024

struct Csc_matrix {
    int columns;
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

__thread_local crts_rply_t dma_rply = 0;
__thread_local unsigned int D_COUNT = 0;
// __thread_local double 

void CMprocess(void *para) {
    CmPara cmpara;
    //split rules
    int *P;
    //start column and end column number
    int start, end;

    //get CPE id
    int tid = CRTS_tid;

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


    //init x and col_ptr in ldm
    double * x = (double *)ldm_malloc(Inc * sizeof(double));
    int * col_ptr  = (int *)ldm_malloc((Inc + 1) * sizeof(int));
    double * data;
    int len;
    for(int round = 0; round < rounds; ++round) {
        int round_col = start + round * Inc;
        //load x and col_ptr of this round to ldm
        CRTS_dma_iget(x, cmpara.splited_x + round_col, 
            Inc * sizeof(double), &dma_rply);
        CRTS_dma_iget(col_ptr, cmpara.csc_matrix.splited_col_ptr + round_col, 
            (Inc + 1) * sizeof(int), &dma_rply);
        D_COUNT += 2;
        CRTS_dma_wait_value(&dma_rply, D_COUNT);
        
        //init data in ldm(data dependence:col_ptr)
        //! if error, test pin length
        int start_idx = col_ptr[0];
        len = col_ptr[Inc] - start_idx;
        data = (double *)ldm_malloc(len * sizeof(double));

        //load value of this round to ldm
        CRTS_dma_iget(data, cmpara.csc_matrix.value + start_idx,
                len * sizeof(double), &dma_rply);
        ++D_COUNT;
        CRTS_dma_wait_value(&dma_rply, D_COUNT);

        //multiply processing
        int idx = 0;
        for(int i = 0; i < Inc; ++i) {
            int num = col_ptr[i + 1] - col_ptr[i];
            for(int j = 0; j < num; ++j) {
                data[idx] = data[idx] * x[i];
                ++idx;
            }
        }
        //! try to using double buffer here(r&w on the same time)
        CRTS_dma_iput(cmpara.csc_matrix.value + start_idx, data, 
                len * sizeof(double), &dma_rply);
        ++D_COUNT;
        CRTS_dma_wait_value(&dma_rply, D_COUNT);
        ldm_free(data);
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

    int idx = 0;
    for(int i = 0; i < end - left; ++i) {
        int num = col_ptr[i + 1] - col_ptr[i];
        for(int j = 0; j < num; ++j) {
            data[idx] = data[idx] * x[i];
            ++idx;
        }
    }
    CRTS_dma_iput(cmpara.csc_matrix.value + col_ptr[0],
            data, len * sizeof(double), &dma_rply);
    ++D_COUNT;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);
    ldm_free(data);
}

void RAprocess(void *para) {

}
