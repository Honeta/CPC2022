#include <slave.h>
#include <simd.h>
#include "matrix_def.h"

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

__thread_local crts_rply_t reply_col = 0, reply_data = 0;
__thread_local unsigned int REPLY_COL_COUNT = 0;
__thread_local unsigned int REPLY_DATA_COUNT = 0;


// Algorithm 5
void naive_thread_block(int thread_id, int start, int end, int L, int R, int *row_ptr, int *col_idx, double *mtx_val, double *vec_val, double *mtx_ans, double *mid_ans)
{
    int start1 = 0, end1 = 0, num = 0, thread = 0, i = 0, j = 0;
    double sum = 0.0;
    thread = thread_id << 1;
    end1 = row_ptr[start];
    for (i = L; i < end1; i++)
        sum = sum + mtx_val[i] * vec_val[col_idx[i]];
    mid_ans[thread] = sum;
    for (i = start; i < end; i++)
    {
        end1 = row_ptr[i + 1];
        sum = 0.0;
        for (j = row_ptr[i]; j < end1; j++)
            sum = sum + mtx_val[j] * vec_val[col_idx[j]];
        mtx_ans[i] = sum;
    }
    start1 = row_ptr[end];
    sum = 0.0;
    for (i = start1; i < R; i++)
        sum = sum + mtx_val[i] * vec_val[col_idx[i]];
    mid_ans[thread + 1] = sum;
}


double calc(int start, int end, int *col_idx, double *data, double *vec) {
    
    int num = end - start;
    double ldm_vec[num] __attribute__((aligned(64)));

    for(int i = 0; i < num; i++) {
        ldm_vec[i] = vec[col_idx[start+i]];
    }
    double sum = .0;
    for(int i = 0; i < num; i++) {
        sum = sum + ldm_vec[i] * data[start+i];
    }

    return sum;
}

// Algorithm 6
void thread_block(int thread_id, int start, int end, int L, int R, int *row_ptr, int *col_idx, double *mtx_val, double *vec_val, double *mtx_ans, double *mid_ans)
{
    int start1 = 0, end1 = 0, num = 0, thread = 0, i = 0;
    thread = thread_id << 1;

    for (i = start; i < end; i++)
    {
        start1 = row_ptr[i];
        end1 = row_ptr[i + 1];
        num = end1 - start1;
        mtx_ans[i] = calc(start1, end1, col_idx, mtx_val, vec_val);
    }
}

double calc_dma_simd(int start, int end, int *col_idx, double *data, double *vec)
{
    int num = end - start;
    int pad = 0;
    if(num % 8 != 0) {
        int div = num >> 3;
        pad = (div<<3)+8-num;
    }
    double ldm_vec[num+pad] __attribute__((aligned(8)));

    for (int i = 0; i < num; i++)
    {
        ldm_vec[i] = vec[col_idx[i]];
    }
    for(int i = num; i < num + pad; i++) {
        ldm_vec[i] = 0;
    }
    double sum = .0;
    doublev8 datav8, vecv8, mulv8;
    for (int i = 0; i < num+pad; i+=8)
    {
        simd_load(datav8, data+i);
        simd_load(vecv8, ldm_vec+i);
        mulv8 = simd_vmuld(datav8, vecv8);
        for(int j = 0; j < 8; j++) {
            sum = sum + mulv8[j];
        }
    }
    return sum;
}

double calc_dma(int start, int end, int *col_idx, double *data, double *vec) {
    

    int num = end - start;
    double ldm_vec[num] __attribute__((aligned(64)));

    for (int i = 0; i < num; i++)
    {
        ldm_vec[i] = vec[col_idx[i]];
    }
    double sum = .0;
    for (int i = 0; i < num; i++)
    {
        sum = sum + ldm_vec[i] * data[i];
    }

    return sum;
}


__thread_local int ldm_col_idx[2][400] __attribute__((aligned(64)));
__thread_local double ldm_data[2][400] __attribute__((aligned(64)));

void thread_block_dma(int thread_id, int start, int end, int L, int R, int *row_ptr, int *col_idx, double *mtx_val, double *vec_val, double *mtx_ans, double *mid_ans)
{
    int start1 = 0, end1 = 0, num = 0, i = 0;

    int lines = end - start;
    int cur = 0, next = 1;
    int ldm_row_ptr[lines + 1] __attribute__((aligned(64)));
    double ldm_ans[lines] __attribute__((aligned(64)));

    CRTS_dma_get(ldm_row_ptr, row_ptr + start, (lines + 1) * sizeof(int));

    int next_num = ldm_row_ptr[1] - ldm_row_ptr[0];

    CRTS_dma_iget(&ldm_col_idx[0][0], col_idx + ldm_row_ptr[0], next_num * sizeof(int), &reply_col);
    CRTS_dma_iget(&ldm_data[0][0], mtx_val + ldm_row_ptr[0], next_num * sizeof(double), &reply_data);
    REPLY_COL_COUNT ++;
    REPLY_DATA_COUNT ++;

    for (i = 0; i < lines - 1; i++)
    {
        CRTS_dma_wait_value(&reply_col, REPLY_COL_COUNT);
        CRTS_dma_wait_value(&reply_data, REPLY_DATA_COUNT);
        next_num = ldm_row_ptr[i+2] - ldm_row_ptr[i+1];
        CRTS_dma_iget(&ldm_col_idx[next][0], col_idx + ldm_row_ptr[i+1], next_num * sizeof(int), &reply_col);
        CRTS_dma_iget(&ldm_data[next][0], mtx_val + ldm_row_ptr[i+1], next_num * sizeof(double), &reply_data);
        REPLY_COL_COUNT++;
        REPLY_DATA_COUNT++;
        start1 = ldm_row_ptr[i];
        end1 = ldm_row_ptr[i + 1];
        ldm_ans[i] = calc_dma_simd(start1, end1, &ldm_col_idx[cur][0], &ldm_data[cur][0], vec_val);
        cur = next;
        next = next ^ 1;
    }

    CRTS_dma_wait_value(&reply_col, REPLY_COL_COUNT);
    CRTS_dma_wait_value(&reply_data, REPLY_DATA_COUNT);
    ldm_ans[lines-1] = calc_dma_simd(ldm_row_ptr[lines-1], ldm_row_ptr[lines], &ldm_col_idx[cur][0], &ldm_data[cur][0], vec_val);
    CRTS_dma_put(mtx_ans + start, ldm_ans, lines * sizeof(double));
}


//Algorithm 10
void slave_func_albus_part1(void* para)
{
    struct AlbusArgs *p = (struct AlbusArgs *)para;
    int tid = CRTS_tid;
    //thread_block(tid, p->start[tid] + 1, p->end[tid], p->block_size[tid], p->block_size[tid+1], p->row_ptr, p->col_idx, p->mtx_val, p->vec_val, p->mtx_ans, p->mid_ans);
    thread_block_dma(tid, p->start[tid] + 1, p->end[tid], p->block_size[tid], p->block_size[tid+1], p->row_ptr, p->col_idx, p->mtx_val, p->vec_val, p->mtx_ans, p->mid_ans);
    return;
}

void slave_func_albus_part2(void* para)
{
    struct AlbusArgs *p = (struct AlbusArgs *)para;
    int tid = CRTS_tid;
    if(tid == 0)
    {
        p->mtx_ans[0] = calc(p->row_ptr[0], p->row_ptr[1], p->col_idx, p->mtx_val, p->vec_val);
        return;
    }
    int sub = tid << 1;
    if(p->end[tid - 1] == p->start[tid]) {
        // p->mtx_ans[p->start[tid]] = p->mid_ans[sub - 1] + p->mid_ans[sub];
        int row = p->start[tid];
        p->mtx_ans[row] = calc(p->row_ptr[row], p->row_ptr[row+1], p->col_idx, p->mtx_val, p->vec_val);
    } else {
        int start = p->start[tid];
        int end = p->end[tid - 1];
        p->mtx_ans[start] = calc(p->row_ptr[start], p->row_ptr[start+1], p->col_idx, p->mtx_val, p->vec_val);
        p->mtx_ans[end] = calc(p->row_ptr[end], p->row_ptr[end+1], p->col_idx, p->mtx_val, p->vec_val);
    }
    return;
}