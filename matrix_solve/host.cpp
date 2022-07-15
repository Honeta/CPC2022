#include <athread.h>
#include "matrix_def.h"
#include "host.h"
#include <math.h>

extern "C" {
    void slave_func(void *);
    void slave_func_cm(void *);
    void slave_func_ra(void *);
}

struct CscMatrix {
    int column;		//列数
    int *column_off;	//非零元在data,col中的索引
    int *rows;		//行标
    double *data;	//数据
    int data_size;
};
int max(int a,int b){return a<=b?b:a;}
int auto_tuner(int NNZ, int N, int L, int *D)
{
    double x = (double)NNZ/N;
    double sum = 0;
    for (int l = 1; l <= L; l++)
    {
        sum = sum + D[l]*(l-x)*(l-x);
    }
    double sigma = sqrt(sum/N);
    int K = L;
    int *DD;
    DD = (int*)malloc(L * sizeof(int));
    DD[0]=0;
    for(int i = 1; i <= L; i++)
    {
        for(int j = 1; j<=L; j++)
        {
            DD[j]=D[j];
        }
        int count = 0;
        for(int j = i+1; j <= L; j++)
        {
            DD[i-1] = DD[i-1] + j/i*D[j];
            if(j%i == 0)
                count = count + (j/i-1)*D[j];
            else
            {
                count = count + (j/i)*D[j];
                DD[j%i-1] = DD[j%i-1] + DD[j-1];
            }
        }
        double xx = NNZ/(N+count);
        double summ = 0;
        for(int l = 1; l <= i; l++)
        {
            summ = summ + DD[l-1]*(l-xx)*(l-xx);
        }
        double sigma2 = sqrt(summ/(N+count));
        if(sigma2 <= sigma/* && enough storage available*/)
            K = i;
    }
    return max(2,K);
}

struct Para_cm{
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
};

void spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    // 实现样例
    // for(int i = 0; i < csr_matrix.rows; i++) {
    //     int start = csr_matrix.row_off[i];
    //     int num = csr_matrix.row_off[i+1] - csr_matrix.row_off[i];
    //     double result = .0;
    //     for(int j = 0; j < num; j++) {
    //         result += x[csr_matrix.cols[start+j]] * csr_matrix.data[start+j];
    //     }
    //     b[i]=result;
    // }
    CRTS_init();
    int NC = CRTS_athread_get_max_threads();
    int NNZ = csr_matrix.row_off[csr_matrix.rows];
    int N = csr_matrix.rows;
    /* csr2csc */
    int column = N;
    int *column_off, *rows;
    double *data;
    column_off = (int*)calloc((N+1),sizeof(int));
    rows = (int*)calloc(NNZ,sizeof(int));
    data = (double*)malloc(NNZ*sizeof(double));
    for(int i = 0; i < NNZ; i++)
    {
        column_off[csr_matrix.cols[i]+1]++;
    }
    for(int i = 0; i < N; i++)
    {
        column_off[i+1] += column_off[i];
    }
    for(int i=0, *ptr=csr_matrix.row_off; i<N; i++, ptr++)
    {
		for (int j=*ptr; j<*(ptr+1); j++){
			int k = csr_matrix.cols[j];
			int l = column_off[k]++;
            rows[l] = i;
			data[l] = csr_matrix.data[j];
		}
	}
    for(int i=N; i>0; i--)
    {
		column_off[i] = column_off[i-1];
	}
	column_off[0] = 0;
    /* split the subA */
    int *D = (int*)calloc((column+1),sizeof(int));
    int L = 0;
    for(int i = 1; i <=column; i++)
    {
        D[column_off[i]-column_off[i-1]]++;
        L = max(L, column_off[i]-column_off[i-1]);
    }
    int K1 = auto_tuner(NNZ, N, L, D);
    int *subA = (int*)malloc((NNZ/K1+1)*sizeof(int));
    double *expand_x = (double*)malloc((NNZ/K1+1)*sizeof(double));
    int ind = 0;
    subA[0]=0;
    expand_x[0] = x[0];
    for(int i = 0; i < N; i++)
    {
        int num = column_off[i+1] - column_off[i];
        int start = column_off[i];
        while(num > K1)
        {
            subA[++ind] = start + K1;
            start += K1;
            num -= K1;
            expand_x[ind] = x[i];
        }
        subA[++ind] = start + num;
        expand_x[ind] = x[i];
    }
    int NNZ_CPE = NNZ / 64;
    int *P = (int*)malloc(65*sizeof(int));
    P[0]=0;
    int CPE_num = 0, tot = 0;
    for(int i = 0; i < ind; i++)
    {
        int num = subA[i+1]-subA[i];
        if(CPE_num >= NNZ_CPE || i == ind -1)
        {
            P[++tot] = i;
            CPE_num = 0;
        }
        else CPE_num += num;
    }
    /* end split */
    double *vals = (double*)calloc(ind, sizeof(double));
    struct Para_cm para;
    para.expand_x = expand_x;
    para.subA = subA;
    para.P = P;
    para.vals = vals;
    para.inc = 1;
    CRTS_athread_spawn(slave_func_cm, &para);
    CRTS_athread_join();
    /////
    CRTS_athread_spawn(slave_func_ra, &para);
    CRTS_athread_join();

    CRTS_athread_halt();
}
