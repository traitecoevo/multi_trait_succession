//*******************************************************************************************/
// HEAD FILE: MATRIX.h
// THIS IS HEAD FILE OF PROGRAM MATRIX.c
//*******************************************************************************************/
# ifndef _MATRIXCPP_H_
# define _MATRIXCPP_H_
# include "utils.h"


/*Note to avoid crashes & program errors, make sure pointers are set to NUll first before calling
these function*/

// 1D, 2D, 3D, AND 4D Int MATRIX MEMORY ALLOCATE
       int *IntV_alloc(int n1);
       int **IntM2D_alloc(int n1, int n2);

// 1D, 2D, 3D, AND 4D DOUBLE MATRIX MEMORY ALLOCATE
       double *Vd_alloc(int n1);
       double *Vd_resize(double *current, int n_old, int n_new);
       double *Vd_expandRow(double *current, int n1, int insert);
       double *Vd_eraseRow(double *current, int n1, int erase);

       double **M2Dd_alloc(int n1, int n2);
       double **M2Dd_resize(double **current, int row_old, int row_new, int col);
       double **M2Dd_eraseRow(double **current, int n1, int n2, int erase);
       double **M2Dd_expandRow(double **current, int n1, int n2, int insert);

// 1D, 2D, 3D, AND 4D FLOAT MATRIX MEMORY ALLOCATE
       float *V_alloc(int n1);
       float **M2D_alloc(int n1, int n2);

// 1D, 2D, 3D, AND 4D FLOAT MATRIX MEMORY ALLOCATE
       int *Vi_alloc(int n1);
       int **M2I_alloc(int n1, int n2);

// MATRIX MEMORY FREE
       void M2DdFree(double **matrix2Dd);
       void M2DFree(float **matrix2D);
       void M2DIFree(int **matrix);

// Allocate memory in arrays from
       double *vector_NR(long nl, long nh);
       void free_vector_NR(double *v, long nl, long nh);
       void nrerror(std::string error_text);
# endif
