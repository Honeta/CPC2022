#ifndef _MATRIX_DEF_H_
#define _MATRIX_DEF_H_

#define MATRIX_NUM  9
const int MESH_NUM[MATRIX_NUM]  = {60000,  60000,  60000, 300000, 300000, 300000, 1200000,  1200000,  1200000 };  	//总网格数
const int ROW_NUM[MATRIX_NUM]	= {10000,  10000,  10000, 50000,  50000,  50000,  200000,   200000,   200000  }; 	//矩阵规模/矩阵行数
const int DATA_NUM[MATRIX_NUM]  = {10000,  30000,  60000, 50000,  150000, 300000, 200000,   600000,   1200000 };	//上三角/下三角非零元个数


struct CsrMatrix {
    int rows;		// 行数
    int *row_off;	// 每行第一个非0元在所有非0元中的索引
    int *cols;		// 第n个非0元素所在的column (matrix * vector, 仅需data所在columns)
    double *data;	// 数据
    int data_size;
};

struct LduMatrix {
	double *upper;
	double *lower;
	double *diag;
	int *uPtr;
	int *lPtr;
	int faces;
	int cells;
};

#include <time.h>
#include <stdio.h>
extern int my_rank;
#define INFO(M, ...) {if(my_rank == 0) { \
                        time_t t; \
                        struct tm *tmif; \
                        t = time(NULL); \
                        tmif = localtime(&t); \
                        printf("[%d-%2d-%2d] [%2d:%2d:%2d] [INFO] " M "",  \
                                tmif->tm_year + 1900, \
                                tmif->tm_mon + 1, \
                                tmif->tm_mday, \
                                tmif->tm_hour, \
                                tmif->tm_min, \
                                tmif->tm_sec, ##__VA_ARGS__); \
                        } \
                    }

#endif