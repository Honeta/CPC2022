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
    int s1 = start1 + 1, s2 = start1 + 2, s3 = start1 + 3;
    doublev8 mtx_3;
    doublev8 mtx_4;
    doublev8 vec_3;
    doublev8 vec_4;
    doublev8 ans_3;
    //Algorithm 7
    if(num >= 4)
    {
        ans = 0.0;
        if(num == 4)
        {
            simd_loadu(mtx_3, mtx_val + start1);
            simd_loadu(mtx_4, mtx_val + s2);
            vec_3 = simd_set_doublev8(vec_val[col_idx[s1]], vec_val[col_idx[start1]],0,0,0,0,0,0);
            vec_4 = simd_set_doublev8(vec_val[col_idx[s3]], vec_val[col_idx[s2]],0,0,0,0,0,0);
            ans_3 = simd_vmuld(mtx_3, vec_3);
            ans_3 = simd_vmad(mtx_4, vec_4, ans_3);
            ans = ans_3[0] + ans_3[1];
        }
        else
        {
            int t = (num>>2)<<2;
            int j_end = start1 + t, num_1 = num & 3;
            int start2 = j_end;
            simd_loadu(mtx_3, mtx_val + start1);
            vec_3 = simd_set_doublev8(vec_val[col_idx[s3]], vec_val[col_idx[s2]], vec_val[col_idx[s1]], vec_val[col_idx[start1]],0,0,0,0);
            doublev8 mtx_ans_1 = simd_vmuld(mtx_3, vec_3);
            start1 = start1 + 4;
            for(int j = start1; j < j_end; j = j + 4)
            {
                s1 = j + 1; s2 = j + 2; s3 = j + 3;
                simd_loadu(mtx_3, mtx_val + j);
                vec_3 = simd_set_doublev8(vec_val[col_idx[s3]], vec_val[col_idx[s2]], vec_val[col_idx[s1]], vec_val[col_idx[j]],0,0,0,0);
                mtx_ans_1 = simd_vmad(mtx_3, vec_3, mtx_ans_1);
            }
            s1 = start2 + 1; s2 = start2 + 2;
            if(num_1 == 0)
            {
                mtx_ans_1 = simd_vaddd(mtx_ans_1, mtx_ans_1);
                ans = ans_3[0] + ans_3[1];
            }
            else if(num_1 == 1)
            {
                mtx_ans_1 = simd_vaddd(mtx_ans_1, mtx_ans_1);
                ans = mtx_ans_1[0]+mtx_ans_1[2]+(mtx_val[start2]*vec_val[col_idx[start2]]);
            }
            else if(num_1 == 2)
            {
                mtx_ans_1 = simd_vaddd(mtx_ans_1, mtx_ans_1);
                ans = mtx_ans_1[0]+mtx_ans_1[2]+(mtx_val[start2]*vec_val[col_idx[start2]]+mtx_val[s1]*vec_val[col_idx[s1]]);
            }
            else
            {
                simd_loadu(mtx_3, mtx_val + start2);
                vec_3 = simd_set_doublev8(0, vec_val[col_idx[s2]], vec_val[col_idx[s1]], vec_val[col_idx[start2]],0,0,0,0);
                mtx_ans_1 = simd_vmad(mtx_3, vec_3, mtx_ans_1);
                mtx_ans_1 = simd_vaddd(mtx_ans_1, mtx_ans_1);
                ans = mtx_ans_1[0]+mtx_ans_1[2];
            }
        }
        return ans;
    }
    //Algorithm 8
    else
    {
        if (num == 0)
            ans = 0.0;
        else
        {
            if(num == 1)
                ans = mtx_val[start1] * vec_val[col_idx[start1]];
            if(num == 2)
                ans = mtx_val[start1]*vec_val[col_idx[start1]]+mtx_val[s1]*vec_val[col_idx[s1]];
            else
            {
                simd_loadu(mtx_3, mtx_val + start1);
                simd_loadu(mtx_4, mtx_val + s2);
                vec_3 = simd_set_doublev8(vec_val[col_idx[s1]],vec_val[col_idx[start1]],0,0,0,0,0,0);
                vec_4 = simd_set_doublev8(0,vec_val[col_idx[s2]],0,0,0,0,0,0);
                ans_3 = simd_vmuld(mtx_3, vec_3);
                ans_3 = simd_vmad(mtx_4, vec_4, ans_3);
                ans = ans_3[0] + ans_3[1];
            }
        }
        return ans;
    }
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
    thread_block_1(tid, p->start[tid] + 1, p->end[tid], p->block_size[tid], p->block_size[tid+1], p->row_ptr, p->col_idx, p->mtx_val, p->vec_val, p->mtx_ans, p->mid_ans);
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