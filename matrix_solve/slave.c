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

double calc(int start, int end, int *col_idx, double *data, double *vec) {
    int num = end - start;
    double ldm_vec[num];
    for(int i = 0; i < num; i++) {
        ldm_vec[i] = vec[col_idx[i+start]];
    }
    double sum = .0;
    for(int i = 0; i < num; i++) {
        sum = sum + ldm_vec[i] * data[i+start];
    }

    return sum;
}

//Algorithm 6
void thread_block(int thread_id, int start, int end, int L, int R, int *row_ptr, int *col_idx, double *mtx_val, double *vec_val, double *mtx_ans, double *mid_ans)
{
    int start1 = 0, end1 = 0, num = 0, thread = 0, i = 0;
    start1 = L;
    end1 = row_ptr[start];
    num = end1 - start1;
    thread = thread_id << 1;
    mid_ans[thread] = calc(start1, end1, col_idx, mtx_val, vec_val);
    for(i = start; i < end; i++)
    {
        start1 = row_ptr[i];
        end1 = row_ptr[i + 1];
        num = end1 - start1;
        mtx_ans[i] = calc(start1, end1, col_idx, mtx_val, vec_val);
    }
    start1 = row_ptr[end];
    end1 = R;
    num = end1 - start1;
    mid_ans[thread + 1] = calc(start1, end1, col_idx, mtx_val, vec_val);
    return;
}
//Algorithm 5
void naive_thread_block(int thread_id, int start, int end, int L, int R, int *row_ptr, int *col_idx, double *mtx_val, double *vec_val, double *mtx_ans, double *mid_ans)
{
    int start1 = 0, end1 = 0, num = 0, thread = 0, i = 0, j = 0;
    double sum = 0.0;
    thread = thread_id << 1;
    end1 = row_ptr[start];
    for(i = L; i < end1; i++)
        sum = sum + mtx_val[i]*vec_val[col_idx[i]];
    mid_ans[thread] = sum;
    for(i = start; i < end; i++)
    {
        end1 = row_ptr[i + 1];
        sum = 0.0;
        for(j = row_ptr[i]; j < end1; j++)
            sum = sum + mtx_val[j]*vec_val[col_idx[j]];
        mtx_ans[i] = sum;
    }
    start1 = row_ptr[end];
    sum = 0.0;
    for(i = start1; i < R; i++)
        sum = sum + mtx_val[i]*vec_val[col_idx[i]];
    mid_ans[thread + 1] = sum;
}

//Algorithm 10
void slave_func_albus_part1(void* para)
{
    struct AlbusArgs *p = (struct AlbusArgs *)para;
    int tid = CRTS_tid;
    //thread_block(tid, p->start[tid] + 1, p->end[tid], p->block_size[tid], p->block_size[tid+1], p->row_ptr, p->col_idx, p->mtx_val, p->vec_val, p->mtx_ans, p->mid_ans);
    thread_block(tid, p->start[tid] + 1, p->end[tid], p->block_size[tid], p->block_size[tid+1], p->row_ptr, p->col_idx, p->mtx_val, p->vec_val, p->mtx_ans, p->mid_ans);
    return;
}

void slave_func_albus_part2(void* para)
{
    struct AlbusArgs *p = (struct AlbusArgs *)para;
    int tid = CRTS_tid;
    if(tid == 0)
    {
        p->mtx_ans[0] = p->mid_ans[0];
        return;
    }
    int sub = tid << 1;
    if(p->end[tid - 1] == p->start[tid]) {
        // p->mtx_ans[p->start[tid]] = p->mid_ans[sub - 1] + p->mid_ans[sub];
        int row = p->start[tid];
        int start = p->row_ptr[row];
        int end = p->row_ptr[row+1];
        double sum = .0;
        for(int i = start; i < end; i++) {
            sum = sum + p->vec_val[p->col_idx[i]] * p->mtx_val[i];
        }
        p->mtx_ans[row] = sum;
    }
    else
    {
        p->mtx_ans[p->start[tid]] = p->mid_ans[sub];
        p->mtx_ans[p->end[tid - 1]] = p->mid_ans[sub - 1];
    }
    return;
}