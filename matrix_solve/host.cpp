#include <athread.h>
#include "matrix_def.h"
#include "host.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ldm_size 256*1024
#define MAX_ITERATIONS 200

struct Csc_matrix {
    int n;              //number of columns
    int *col_ptr;       //index of the first element in a column
    int *col_emt;       //number of element in every column
    int *rows;          //row number of every non-0 element
    double *value;      //value of every non-0 element
    int data_size;

    int inc_cnt;
    int *splited_col_ptr;
};

struct CMPara {
    struct Csc_matrix csc_matrix;
    int *P;
    double *splited_x;
    int Inc;
};

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

struct RAPara {
    struct Result_matrix result_matrix;
    int *P;
    double *splited_b;
    int Inc;
};

//threshold
int K_cm = 0, K_ra = 0;
//CPE level partition pointer
//elements of P_cm is column index, elements of P_ra is row index
int *P_cm, *P_ra;
//mapping from x to splited x, mapping from splited b to b
int *splited_x_map, *splited_b_map;
//expanded x, b before reduce
double *splited_x, *splited_b;
//store static csc value for every iterate
double *static_value;
//global variable of iter round
int iter_rounds = 0;
Csc_matrix csc_matrix;
Result_matrix result_matrix;
CMPara para_cm;
RAPara para_ra;


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
    //number of element in every column
    for (int i = 0; i < data_size; ++i) {
        ++c_col_ptr[r_cols[i]];
    }
    //store the number of element of every column to col_emt
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
    csc_matrix.n         = n;
    csc_matrix.data_size = data_size;
    csc_matrix.col_ptr   = c_col_ptr;
    csc_matrix.col_emt   = col_emt;
    csc_matrix.rows      = c_rows;
    csc_matrix.value     = (double *)malloc(data_size * sizeof(double));
}

//convert from csc format to csr format
void csc2csr(const CsrMatrix &csr_matrix) {
    int data_size   = csc_matrix.data_size;
    int *rows       = csc_matrix.rows;
    double *c_data  = csc_matrix.value;
    int *row_ptr    = result_matrix.row_ptr;
    double *r_data  = result_matrix.value;
    for (int i = 0; i < data_size; ++i) {
        int row = rows[i];
        r_data[row_ptr[row]] = c_data[i];
        ++row_ptr[row];
    }
    memcpy(result_matrix.row_ptr, csr_matrix.row_off, csr_matrix.rows * sizeof(int));
}

//init columns multiply result matrix in csr format
void result_init(const CsrMatrix &csr_matrix) {
    result_matrix.n       = csr_matrix.rows;
    result_matrix.row_ptr = (int *)malloc(csr_matrix.rows * sizeof(int));
    result_matrix.columns = (int *)malloc(csr_matrix.data_size * sizeof(int));
    result_matrix.value   = (double *)malloc(csr_matrix.data_size * sizeof(double));
    result_matrix.data_size = csr_matrix.data_size;
    memcpy(result_matrix.row_ptr, csr_matrix.row_off, csr_matrix.rows * sizeof(int));
    memcpy(result_matrix.columns, csr_matrix.cols, csr_matrix.data_size * sizeof(int));
}

//calculate new threshold K_cm of columns vector for CM
void auto_tuner_cm() {
    int n = csc_matrix.n;
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
    K_cm = L;
    for (int i = 2; i < L; ++i) {
        int * DD = (int *)malloc((i+1) * sizeof(int));
        memcpy(DD, D, (i+1) * sizeof(int));
        //split long vec into short vec in statistic
        for (int j = i + 1; j <= L; ++j) {
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
            K_cm = i;
            csc_matrix.inc_cnt = inc_cnt;
        }
        free(DD);
    }
    free(D);
}

//split vec into vec whose length more than K_cm and record P
void vec_split_cm() {
    int n = csc_matrix.n;
    int inc_cnt = csc_matrix.inc_cnt;
    int *col_ptr = csc_matrix.col_ptr;
    int *col_emt = csc_matrix.col_emt;
    int *splited_col_ptr = 
        (int *)malloc((n + inc_cnt) * sizeof(int));
    splited_x_map = 
        (int *)malloc((n + inc_cnt) * sizeof(int));
    splited_x = (double *)malloc((n + inc_cnt) * sizeof(double));
    P_cm = (int *)malloc(65 * sizeof(int));

    int idx_ptr = 0;
    for (int i = 0; i < n; ++i) {
        int count = col_emt[i] / K_cm + 1;
        splited_col_ptr[idx_ptr] = col_ptr[i];
        for (int j = 0; j < count; ++j) {
            if(j != 0) {
                splited_col_ptr[idx_ptr] = splited_col_ptr[idx_ptr - 1] + K_cm;
            }
            if(idx_ptr % ((n + inc_cnt) / 64) == 0) {
                P_cm[idx_ptr / ((n + inc_cnt) / 64)] = idx_ptr;
            }
            splited_x_map[idx_ptr] = i;
            ++idx_ptr;
        }
    }
    P_cm[65] = n + inc_cnt;
    csc_matrix.splited_col_ptr = splited_col_ptr;
}

//calculate new threshold K_ra of columns vector for RA
void auto_tuner_ra(const CsrMatrix &csr_matrix) {
    int n = csr_matrix.rows;
    int data_size = csr_matrix.data_size;
    int * row_ptr = csr_matrix.row_off;
    int * row_emt = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n - 1; ++i) {
        row_emt[i] = row_ptr[i + 1] - row_ptr[i];
    }
    row_emt[n - 1] = data_size - row_ptr[n - 1];
    result_matrix.row_emt = row_emt;
    //calculate the maximum number of elements in a column
    int L = 0;
    for (int i = 0; i < n; ++i) {
        if(row_emt[i] > L) {
            L = row_emt[i];
        }
    }

    //D[i] : the number of columns which include i element
    //calculate D 
    int * D = (int *)malloc(L * sizeof(int));
    for (int i = 0; i < n; ++i) {
        ++D[row_emt[i]];
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
    K_ra = L;
    for (int i = 2; i < L; ++i) {
        int * DD = (int *)malloc((i+1) * sizeof(int));
        memcpy(DD, D, (i+1) * sizeof(int));
        //split long vec into short vec in statistic
        for (int j = i + 1; j <= L; ++j) {
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
            K_ra = i;
            result_matrix.inc_cnt = inc_cnt;
        }
        free(DD);
    }
    free(D);
}

//split vec into vec whose length more than K_ra and record P
void vec_split_ra() {
    int n = result_matrix.n;
    int inc_cnt = result_matrix.inc_cnt;
    int *row_ptr = result_matrix.row_ptr;
    int *row_emt = result_matrix.row_emt;
    int *splited_row_ptr = 
        (int *)malloc((n + inc_cnt) * sizeof(int));
    splited_b_map = 
        (int *)malloc((n + inc_cnt) * sizeof(int));
    splited_b = (double *)malloc((n + inc_cnt) * sizeof(double));
    P_ra = (int *)malloc(65 * sizeof(int));
    int idx_ptr = 0;
    for (int i = 0; i < n; ++i) {
        int count = row_emt[i] / K_ra + 1;
        splited_row_ptr[idx_ptr] = row_ptr[i];
        for (int j = 0; j < count; ++j) {
            if(j != 0) {
                splited_row_ptr[idx_ptr] = splited_row_ptr[idx_ptr - 1] + K_ra;
            }
            if(idx_ptr %((n + inc_cnt) / 64) == 0) {
                P_ra[idx_ptr / ((n + inc_cnt) / 64)] = idx_ptr;
            }
            splited_b_map[idx_ptr] = i;
            ++idx_ptr;
        }
    }
    P_ra[65] = n + inc_cnt;
    result_matrix.splited_row_ptr = splited_row_ptr;
}

//expand vector x to fit in the threshold K
void x_expand(double *x) {
    int length = csc_matrix.n + csc_matrix.inc_cnt;
    for (int i = 0; i < length; ++i) {
        splited_x[i] = x[splited_x_map[i]];
    }
}

void b_reduce(double *b) {
    memset(b, 0, result_matrix.n * sizeof(double));
    int len = result_matrix.n + result_matrix.inc_cnt;
    for (int i = 0; i < len; ++i) {
        b[splited_b_map[i]] += splited_b[i];
    }
}

//release malloced space per MAX_ITERATIONS iterate rounds
void release_space() {
    free(csc_matrix.col_emt);
    free(csc_matrix.col_ptr);
    free(csc_matrix.splited_col_ptr);
    free(csc_matrix.rows);
    free(csc_matrix.value);
    free(result_matrix.row_emt);
    free(result_matrix.splited_row_ptr);
    free(static_value);
    free(splited_x);
    free(splited_x_map);
    free(P_cm);
    free(P_ra);
}

extern "C" {
    void CMprocess(void *para);
    void RAprocess(void *para);
}

//parallel spmv kernel
void spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    if(iter_rounds == 0) {
        CRTS_init();
    }
    //update per MAX_ITERATIONS rounds 
    if(iter_rounds % MAX_ITERATIONS == 0) {
        csr2csc(csr_matrix);
        result_init(csr_matrix);
        auto_tuner_cm();
        vec_split_cm();
        //get K_ra
        auto_tuner_ra(csr_matrix);
        //get P_ra
        vec_split_ra();
        //set columns multiply para
        memcpy(&(para_cm.csc_matrix), &csc_matrix, sizeof(struct Csc_matrix));
        para_cm.splited_x = splited_x;
        para_cm.P = P_cm;
        para_cm.Inc = ldm_size / (8 * K_cm + 12);
        //set rows addtion para
        memcpy(&(para_ra.result_matrix), &result_matrix, sizeof(struct Result_matrix));
        para_ra.splited_b = splited_b;
        para_ra.P = P_ra;
        para_ra.Inc = ldm_size / (8 * K_ra + 12);
    }
    ++iter_rounds;

    //update splited_x every iteration by splited_x_map
    //renew csc_matrix.value every iteration(value also store CM result)
    x_expand(x);
    memcpy(csc_matrix.value, static_value, csc_matrix.data_size * sizeof(double));
    CRTS_athread_spawn_noflush(CMprocess, &para_cm);
    CRTS_athread_join();

    //only need to reposition value, reuse row_ptr and columns of csr_matrix
    csc2csr(csr_matrix);


    CRTS_athread_spawn_noflush(RAprocess, &para_ra);
    CRTS_athread_join();
    b_reduce(b);

    //release space when matrix has been changed
    if(iter_rounds % MAX_ITERATIONS == 0) {
        release_space();
    }
    //release CPE resource after all the calculation has finished
    if(iter_rounds == 9 * MAX_ITERATIONS) {
        CRTS_athread_halt();
    }
}

void naive_spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    // sample
    for (int i = 0; i < csr_matrix.rows; i++) {
        //get the index of the first element in the row
        int start = csr_matrix.row_off[i];
        //get the number of non-zero element in the row
        int num = csr_matrix.row_off[i+1] - csr_matrix.row_off[i];
        double result = .0;
        for (int j = 0; j < num; j++) {
            result += x[csr_matrix.cols[start+j]] * csr_matrix.data[start+j];
        }
        b[i]=result;
    }
}
