#include <athread.h>
#include "matrix_def.h"
#include "host.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ldm_size 256*1024
#define MAX_ITERATIONS 200

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

struct Result_matrix {
    int *row_ptr;
    int *columns;
    double *value;
};

struct CMPara {
    struct Csc_matrix csc_matrix;
    int *P;
    double *splited_x;
    int Inc;
};

//threshold
int K = 0;
//CPE level partition pointer, element is column index 
int *P;
//mapping from x to splited x
//update every 200 rounds
int *splited_x_map;
//expanded x
double *splited_x;
double *static_value, *c_value;
//global variable of iter round
int iter_rounds = 0;
Csc_matrix csc_matrix;
Result_matrix result_matrix;


//convert from csr format to csc format
void csr2csc(const CsrMatrix &csr_matrix) {
    int data_size    = csr_matrix.data_size;
    int n            = csr_matrix.rows;
    int * r_row_ptr  = csr_matrix.row_off;
    int * r_cols     = csr_matrix.cols;
    double * r_value = csr_matrix.data;
    int * c_col_ptr  = (int *)malloc(n * sizeof(int));
    int * col_emt    = (int *)malloc(n * sizeof(int));
    int * c_rows     = (int *)malloc(data_size * sizeof(int));
    static_value = (double *)malloc(data_size * sizeof(double));
    c_value = (double *)malloc(data_size * sizeof(double));
    //number of element in every column
    for (int i = 0; i < data_size; ++i) {
        ++c_col_ptr[r_cols[i]];
    }
    memcpy(col_emt, c_col_ptr, n * sizeof(int));
    //index of the first element in every column
    int sum = 0;
    int tmp = 0;
    for (int col = 0; col < n; ++col) {
        tmp = c_col_ptr[col];
        c_col_ptr[col] = sum;
        sum += tmp;
    }
    //load data into csc value and row array
    for (int row = 0; row < n; ++row) {
        int end = (row == n - 1 ? data_size : r_row_ptr[row + 1]);
        for (int i = r_row_ptr[row]; i < end; ++i) {
            int col = r_cols[i];
            int index = c_col_ptr[col];
            static_value[index] = r_value[i];
            c_rows[index] = row;
            ++c_col_ptr[col];
        }
    }
    //recover c_col_ptr array
    int first = 0;
    for (int col = 0; col < n; ++col) {
        tmp = c_col_ptr[col];
        c_col_ptr[col] = tmp;
        first = tmp;
    }
    //set struct
    csc_matrix.columns   = n;
    csc_matrix.data_size = data_size;
    csc_matrix.col_ptr   = c_col_ptr;
    csc_matrix.col_emt   = col_emt;
    csc_matrix.rows      = c_rows;
    csc_matrix.value     = c_value;
}

//convert from csc format to csr format
void csc2csr(const CsrMatrix &csr_matrix) {
    int data_size   = csc_matrix.data_size;
    int *rows       = csc_matrix.rows;
    double *c_data  = csc_matrix.value;
    int *row_ptr    = result_matrix.row_ptr;
    double *r_data  = result_matrix.value;
    for(int i = 0; i < data_size; ++i) {
        int row = rows[i];
        r_data[row_ptr[row]] = c_data[i];
        ++row_ptr[row];
    }
    memcpy(row_ptr, csr_matrix.row_off, csr_matrix.rows * sizeof(int));
}

//calculate new threshold K
void auto_tuner() {
    int n = csc_matrix.columns;
    int data_size = csc_matrix.data_size;
    int * col_emt = csc_matrix.col_emt;
    //calculate the maximum number of elements in a column
    int L = 0;
    for (int i = 0; i < n; ++i) {
        if(col_emt[i] > L) {
            L = col_emt[i];
        }
    }

    //D[i] : the number of columns which include i element
    //calculate D 
    int * D = (int *)malloc(L * sizeof(int));
    for (int i = 0; i < n; ++i) {
        ++D[col_emt[i]];
    }

    //init std
    double sum = 0;
    double std = 0;
    double std_tmp = 0;
    double mean = data_size / n;
    for (int i = 0; i < L; ++i) {
        sum += D[i] * (i - mean) * (i - mean);
    }
    std = sqrt(sum / n);

    //optmize K to minimum std, K ranging from 2 to L
    int inc_cnt = 0;
    K = L;
    for (int i = 2; i < L; ++i) {
        int * DD = (int *)malloc((i+1) * sizeof(int));
        memcpy(DD, D, (i+1) * sizeof(int));
        //split long vec into short vec in statistic
        for(int j = i + 1; j <= L; ++j) {
            //vector which length is threshold
            DD[i] += (j / i) * D[j];
            //no remainder vector
            if(j % i == 0){
                inc_cnt += D[j] * (j / i - 1);
            }
            //processing remainder vectors
            else{
                DD[j % i] += D[j];
                inc_cnt += D[j] * (j / i);
            }
        }

        //calculate std
        mean = data_size / (n + inc_cnt);
        sum = 0;
        for (int j = 0; j <= i; ++j) {
            sum += DD[j] * (j - mean) * (j - mean);
        }
        std_tmp = sqrt(sum / (n + inc_cnt));
        //update K and std
        if(std_tmp < std) {
            std = std_tmp;
            K = i;
            csc_matrix.inc_cnt = inc_cnt;
        }
        free(DD);
    }
    free(D);
}

//init columns multiply result matrix in csr format
void result_init(const CsrMatrix &csr_matrix) {
    result_matrix.row_ptr = (int *)malloc(csr_matrix.rows * sizeof(int));
    result_matrix.columns = (int *)malloc(csr_matrix.data_size * sizeof(int));
    result_matrix.value   = (double *)malloc(csr_matrix.data_size * sizeof(double));
    memcpy(result_matrix.row_ptr, csr_matrix.row_off, csr_matrix.rows * sizeof(int));
    memcpy(result_matrix.columns, csr_matrix.cols, csr_matrix.data_size * sizeof(int));
}

//split vec into vec whose length less than K
//and record CPE level partition pointer array P
//four-ways partition
void vec_split() {
    int n = csc_matrix.columns;
    int inc_cnt = csc_matrix.inc_cnt;
    int *col_ptr = csc_matrix.col_ptr;
    int *col_cmt = csc_matrix.col_emt;
    int *splited_col_ptr = 
        (int *)malloc((n + inc_cnt) * sizeof(int));
    splited_x_map = 
        (int *)malloc((n + inc_cnt) * sizeof(int));
    P = (int *)malloc(65 * sizeof(int));
    int idx_ptr = 0;
    for (int i = 0; i < n; ++i) {
        int count = col_cmt[i] / K + 1;
        splited_col_ptr[idx_ptr] = col_ptr[i];
        for(int j = 0; j < count; ++j) {
            if(j != 0) {
                splited_col_ptr[idx_ptr] = splited_col_ptr[idx_ptr - 1] + K;
            }
            if(idx_ptr % ((n + inc_cnt) / 64) == 0) {
                P[idx_ptr / ((n + inc_cnt) / 64)] = idx_ptr;
            }
            P[65] = n + inc_cnt;
            splited_x_map[idx_ptr] = i;
            ++idx_ptr;
        }
    }
    csc_matrix.splited_col_ptr = splited_col_ptr;
}

//expand vector x to fit in the threshold K
void x_expand(double *x) {
    int length = csc_matrix.columns + csc_matrix.inc_cnt;
    splited_x = (double *)malloc(length * sizeof(double));
    for(int i = 0; i < length; ++i) {
        splited_x[i] = x[splited_x_map[i]];
    }
}

//release malloced space per MAX_ITERATIONS iterate rounds
void release_space() {
    free(csc_matrix.col_emt);
    free(csc_matrix.col_ptr);
    free(csc_matrix.splited_col_ptr);
    free(csc_matrix.rows);
    free(c_value);
    free(static_value);
    free(splited_x);
    free(splited_x_map);
    free(P);
}

extern "C" {
    void CMprocess(void *CMpara);
}

//parallel spmv kernel
void spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    if(iter_rounds == 0) {
        CRTS_init();
    }
    //update per MAX_ITERATIONS rounds 
    if(iter_rounds % MAX_ITERATIONS == 0) {
        csr2csc(csr_matrix);
        auto_tuner();
        vec_split();
        result_init(csr_matrix);
    }
    ++iter_rounds;
    //update splited_x every iteration(x updated every iteration)
    x_expand(x);
    //renew c_value every iteration(c_value also store CM result)
    memcpy(c_value, static_value, csc_matrix.data_size * sizeof(double));
    

    //set columns multiply para
    CMPara para;
    para.csc_matrix = csc_matrix;
    para.splited_x = splited_x;
    para.P = P;
    para.Inc = ldm_size / (8 * K + 12);

    //Columns Multiply process
    CRTS_athread_spawn_noflush(CMprocess, &para);
    CRTS_athread_join();

    //convert csc format to csr format to accelerate RA
    csc2csr(csr_matrix);
    //set rows addtion para



    //release space when matrix has been changed
    if(iter_rounds % MAX_ITERATIONS == 0) {
        release_space();
    }
    //release CPE resource after all the calculation has finished
    if(iter_rounds == 9 * MAX_ITERATIONS) {
        CRTS_athread_halt();
    }
}

// x : vector
// b : result
void native_spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    // sample
    for(int i = 0; i < csr_matrix.rows; i++) {
        //get the index of the first element in the row
        int start = csr_matrix.row_off[i];
        //get the number of non-zero element in the row
        int num = csr_matrix.row_off[i+1] - csr_matrix.row_off[i];
        double result = .0;
        for(int j = 0; j < num; j++) {
            result += x[csr_matrix.cols[start+j]] * csr_matrix.data[start+j];
        }
        b[i]=result;
    }
}
