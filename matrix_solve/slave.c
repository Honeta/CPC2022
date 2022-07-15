#include <slave.h>
#include "matrix_def.h"
#include <stdio.h>
#include <athread.h>
#include <crts.h>

typedef struct {
    int inc;
    int *P;
    int *subA;
    double *expand_x;
    double *vals;
    double *data;
    int *column_off;
    int *rows;
    int NNZ;
    int N;
}Para_cm;

void slave_func(void* para)
{
    printf("hello\n");
    return;
}
void slave_func_cm(void* para)
{
    crts_rply_t dma_reply = 0;
    unsigned int D_COUNT = 0;
    int my_id = CRTS_rid*8+CRTS_cid;
    
    Para_cm para_s;
    CRTS_dma_iget(&para_s, para, sizeof(Para_cm),&dma_reply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_reply, D_COUNT);
    // __thread_local int P[MAXA] ;
    // __thread_local double vals[MAXA] __attribute__((aligned(64)));
    // __thread_local double seg_x[MAXA] __attribute__((aligned(64)));
    // __thread_local int subA[MAXA];
     int start = para_s.P[my_id];
     int end = para_s.P[my_id+1];
     int i = start;
     while(i < end)
    {
        int num = para_s.subA[i+1] - para_s.subA[i];
        int start2 = para_s.subA[i];
        for(int j = 0; j < num; j++)
        para_s.rows[start2 + j] = para_s.rows[start2 + j]*para_s.expand_x[i+1];
        i = i + para_s.inc;
    }
    return;
}
void slave_func_ra(void* para)
{
    return;
}
