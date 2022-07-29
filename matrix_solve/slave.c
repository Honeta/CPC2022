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

//Algorithm 9
double calculation(int start1, int end1, int num, int *row_ptr, int *col_idx, double *mtx_val, double *vec_val)
{
    double ans;
    int s0 = start1, s1 = start1 + 1, s2 = start1 + 2, s3 = start1 + 3, s4 = start1 + 4, s5 = start1 + 5, s6 = start1 + 6, s7 = start1 + 7;
    doublev8 mtx_8;
    doublev8 vec_8;
    doublev8 ans_8;
    double store_ans[8];//__attribute__ ((aligned (64)));
    //Algorithm 7
    if(num >= 8)
    {
        ans = 0.0;
        if(num == 8)
        {
            simd_loadu(mtx_8, mtx_val + s0);
            vec_8 =
                simd_set_doublev8(vec_val[col_idx[s0]], vec_val[col_idx[s1]],
                                  vec_val[col_idx[s2]], vec_val[col_idx[s3]],
                                  vec_val[col_idx[s4]], vec_val[col_idx[s5]],
                                  vec_val[col_idx[s6]], vec_val[col_idx[s7]]);
            ans_8 = simd_vmuld(mtx_8, vec_8);
            ans = simd_reduc_plusd(ans_8);
        }
        else
        {
            int t = (num>>3)<<3;
            int j_end = start1 + t, num_1 = num & 7;
            int start2 = j_end;
            simd_loadu(mtx_8, mtx_val + s0);
            vec_8 =
                simd_set_doublev8(vec_val[col_idx[s0]], vec_val[col_idx[s1]],
                                  vec_val[col_idx[s2]], vec_val[col_idx[s3]],
                                  vec_val[col_idx[s4]], vec_val[col_idx[s5]],
                                  vec_val[col_idx[s6]], vec_val[col_idx[s7]]);
            doublev8 mtx_ans_1 = simd_vmuld(mtx_8, vec_8);
            start1 = start1 + 8;
            for(int j = start1; j < j_end; j = j + 8)
            {
                s0 = j, s1 = j + 1; s2 = j + 2; s3 = j + 3, s4 = j + 4, s5 = j + 5, s6 = j + 6, s7 = j + 7;
                simd_loadu(mtx_8, mtx_val + j);
                vec_8 = simd_set_doublev8(
                    vec_val[col_idx[s0]], vec_val[col_idx[s1]],
                    vec_val[col_idx[s2]], vec_val[col_idx[s3]],
                    vec_val[col_idx[s4]], vec_val[col_idx[s5]],
                    vec_val[col_idx[s6]], vec_val[col_idx[s7]]);
                mtx_ans_1 = simd_vmad(mtx_8, vec_8, mtx_ans_1);
            }
            s0 = start2, s1 = start2 + 1; s2 = start2 + 2, s3 = start2 + 3, s4 = start2 + 4, s5 = start2 + 5, s6 = start2 + 6, s7 = start2 + 7;
            ans = simd_reduc_plusd(mtx_ans_1);
            for (int k = 0; k < num_1; k++)
                ans = ans + mtx_val[s0 + k] * vec_val[col_idx[s0 + k]];
        }
        return ans;
    }
    //Algorithm 8
    else
    {
        if (num == 0)
            ans = 0.0;
        if (num == 1)
            ans = mtx_val[s0] * vec_val[col_idx[s0]];
        if (num == 2)
            ans = mtx_val[s0] * vec_val[col_idx[s0]] + mtx_val[s1] * vec_val[col_idx[s1]];
        if (num == 3)
            ans = mtx_val[s0] * vec_val[col_idx[s0]] + mtx_val[s1] * vec_val[col_idx[s1]] + mtx_val[s2] * vec_val[col_idx[s2]];
        if (num == 4)
            ans = mtx_val[s0] * vec_val[col_idx[s0]] + mtx_val[s1] * vec_val[col_idx[s1]] + mtx_val[s2] * vec_val[col_idx[s2]] + mtx_val[s3] * vec_val[col_idx[s3]];
        if (num == 5)
            ans = mtx_val[s0] * vec_val[col_idx[s0]] + mtx_val[s1] * vec_val[col_idx[s1]] + mtx_val[s2] * vec_val[col_idx[s2]] + 
            mtx_val[s3] * vec_val[col_idx[s3]] + mtx_val[s4] * vec_val[col_idx[s4]];
        if (num == 6)
            ans = mtx_val[s0] * vec_val[col_idx[s0]] + mtx_val[s1] * vec_val[col_idx[s1]] + mtx_val[s2] * vec_val[col_idx[s2]] + 
            mtx_val[s3] * vec_val[col_idx[s3]] + mtx_val[s4] * vec_val[col_idx[s4]] + mtx_val[s5] * vec_val[col_idx[s5]];
        if (num == 7)
            ans = mtx_val[s0] * vec_val[col_idx[s0]] + mtx_val[s1] * vec_val[col_idx[s1]] + mtx_val[s2] * vec_val[col_idx[s2]] + 
            mtx_val[s3] * vec_val[col_idx[s3]] + mtx_val[s4] * vec_val[col_idx[s4]] + mtx_val[s5] * vec_val[col_idx[s5]] + mtx_val[s6] * vec_val[col_idx[s6]];      
    }
    return ans;
}

//Algorithm 6
void thread_block(int thread_id, int start, int end, int L, int R, int *row_ptr, int *col_idx, double *mtx_val, double *vec_val, double *mtx_ans, double *mid_ans)
{
    int start1 = 0, end1 = 0, num = 0, thread = 0, i = 0;
    start1 = L;
    end1 = row_ptr[start];
    num = end1 - start1;
    thread = thread_id << 1;
    mid_ans[thread] = calculation(start1, end1, num, row_ptr, col_idx, mtx_val, vec_val);
    for(i = start; i < end; i++)
    {
        start1 = row_ptr[i];
        end1 = row_ptr[i + 1];
        num = end1 - start1;
        mtx_ans[i] = calculation(start1, end1, num, row_ptr, col_idx, mtx_val, vec_val);
    }
    start1 = row_ptr[end];
    end1 = R;
    num = end1 - start1;
    mid_ans[thread + 1] = calculation(start1, end1, num, row_ptr, col_idx, mtx_val, vec_val);
    return;
}
//Algorithm 5
void thread_block_1(int thread_id, int start, int end, int L, int R, int *row_ptr, int *col_idx, double *mtx_val, double *vec_val, double *mtx_ans, double *mid_ans)
{
    int start1 = 0, end1 = 0, num = 0, thread = 0, i = 0, j = 0;
    double sum = 0.0;
    thread = thread_id << 1;
    end1 = row_ptr[start];
    for(i = L; i < end1; i++)
        sum += mtx_val[i]*vec_val[col_idx[i]];
    mid_ans[thread] = sum;
    for(i = start; i < end; i++)
    {
        end1 = row_ptr[i + 1];
        sum = 0.0;
        for(j = row_ptr[i]; j < end1; j++)
            sum += mtx_val[j]*vec_val[col_idx[j]];
        mtx_ans[i] = sum;
    }
    start1 = row_ptr[end];
    sum = 0.0;
    for(i = start1; i < R; i++)
        sum += mtx_val[i]*vec_val[col_idx[i]];
    mid_ans[thread + 1] = sum;
}

//Algorithm 10
void slave_func_albus_part1(void* para)
{
    struct AlbusArgs *p = (struct AlbusArgs *)para;
    int tid = CRTS_tid;
    thread_block(tid, p->start[tid] + 1, p->end[tid], p->block_size[tid], p->block_size[tid+1], p->row_ptr, p->col_idx, p->mtx_val, p->vec_val, p->mtx_ans, p->mid_ans);
    //thread_block_1(tid, p->start[tid] + 1, p->end[tid], p->block_size[tid], p->block_size[tid+1], p->row_ptr, p->col_idx, p->mtx_val, p->vec_val, p->mtx_ans, p->mid_ans);
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
    if (p->end[tid - 1] == p->start[tid])
    {
      // p->mtx_ans[p->start[tid]] = p->mid_ans[sub - 1] + p->mid_ans[sub];
      int row = p->start[tid];
      int start = p->row_ptr[row];
      int end = p->row_ptr[row + 1];
      double sum = .0;
      for (int i = start; i < end; i++) {
        sum = sum + p->vec_val[p->col_idx[i]] * p->mtx_val[i];
      }
      p->mtx_ans[row] = sum;
    } else {
      p->mtx_ans[p->start[tid]] = p->mid_ans[sub];
      p->mtx_ans[p->end[tid - 1]] = p->mid_ans[sub - 1];
    }
    return;
}