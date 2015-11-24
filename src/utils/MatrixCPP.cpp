//*******************************************************************************************/
//
//    PROGRAM NAME: MATRIX.c
//    PURPOSE: ALLOCATE DYNAMIC MEMORY FOR VECTOR, 2D, 3D, AND 4D MATRIX.
//
//    INCLUDE HEAD FILE: MATRIX.h
//
//    SUB-PROGRAMS FOR CALLING IN THIS C FILE ARE AS FOLLOWS:
//
//      I. DOUBLE FLOAT MATRIX MEMORY ALLOCATION:
//         1. double *Vd_alloc(int n1): FLOAT VECTOR
//         2. double **M2Dd_alloc(int n1, int n2): 2 DIMENSIONAL FLOAT MATRIX
//         3. double ***M3Dd_alloc(int n1, int n2, int n3): 3 DIMENSIONAL FLOAT MATRIX
//         4. double ****M4Dd_alloc(int n1, int n2, int n3, int n4): 4 DIMENSIONAL FLOAT MATRIX
//
//     II. FLOAT MATRIX MEMORY ALLOCATION:
// 5. float *V_alloc(int n1): FLOAT VECTOR
// 6. float **M2D_alloc(int n1, int n2): 2 DIMENSIONAL FLOAT MATRIX
// 7. float ***M3D_alloc(int n1, int n2, int n3): 3 DIMENSIONAL FLOAT MATRIX
// 8. float ****M4D_alloc(int n1, int n2, int n3, int n4): 4 DIMENSIONAL FLOAT MATRIX
//
//    III. INTEGER MATRIX MEMORY ALLOCATION:
// 9. int *IntV_alloc(int n1): INTEGER VECTOR
//        10. int **IntM2D_alloc(int n1, int n2): 2 D. INTEGER MATRIX
//        11. int ***IntM3D_alloc(int n1, int n2, int n3): 3 D. INTEGER MATRIX
//        12. int ****IntM4D_alloc(int n1, int n2, int n3, int n4): 4 D. INTEGER MATRIX
//
//        13-1. long *IntVd_alloc(int n1): INTEGER VECTOR
//        13-2. long *IntM2Dd_alloc(int n1, int n2): INTEGER 2-D MATRIX
//
//     IV. FREE MEMORY:
//        14  free( *vector): 1 D. VECTOR MEMORY FREE
//        15. M2DdFree(double **matrix2D, int n1): 2 D. MATRIX MEMORY FREE
//        16. M3DdFree(double ***matrix3D, int n1, int n2): 3D. MATRIX MEMORY FREE
//        17. M4DdFree(double ****matrix4D, int n1, int n2, int n3): 4 D. MATRIX MEMORY FREE
//        18. M2DFree(float **matrix2D, int n1): 2 D. MATRIX MEMORY FREE
//        19. M3DFree(float ***matrix3D, int n1, int n2): 3D. MATRIX MEMORY FREE
//        20. M4DFree(float ****matrix4D, int n1, int n2, int n3): 4 D. MATRIX MEMORY FREE
//        21. IntM2DFree(int **matrix2D, int n1): 2 D. MATRIX MEMORY FREE
//        22. INtM3DFree(int ***matrix3D, int n1, int n2): 3D. MATRIX MEMORY FREE
//        23. IntM4DFree(int ****matrix4D, int n1, int n2, int n3): 4 D. MATRIX MEMORY FREE
//        24. IntM2DdFree(long **matrix2D, int n1): 2 D. MATRIX MEMORY FREE
//*******************************************************************************************/

# include "MatrixCPP.h"
#include <string.h> // Needed for NULL pointer

#include <iostream>
#include <fstream>
# include <stdlib.h>

//*******************************************************************************************/
//     DOUBLE FLOAT MATRIX ALLOCATION:
// 1. DOUBLE VECTOR
// 2. DOUBLE MATRIX2D
// 3. DOUBLE MATRIX3D
// 4. DOUBLE MATRIX4D
//*******************************************************************************************/
/* 1. DOUBLE VECTOR */
double *Vd_alloc(int n1) {
	double *V = new double[n1];
	return V;
}

double *Vd_resize(double *current, int n_old, int n_new) {
	double *a;
	if(n_new>n_old)
		a=Vd_expandRow(current, n_old, n_new- n_old);
	else
	{
		int i=0;
		a =Vd_alloc(n_new);
				// Copy data before erase row
		while(i<n_new)
		{
			a[i] = current[i];
			i++;
		}
		delete [] current;
	}
	return a;
}


/*Removes row form existing vector*/
double* Vd_eraseRow(double *current, int n1, int erase) {
	int i=0;
	double *a =Vd_alloc(n1-1);
	// Copy data before erase row
	while(i<erase)
	{
		a[i] = current[i];
		i++;
	}
	// Skip row = erase
	// Copy data after erased row
	i = erase+1;
	while(i<n1)
	{
		a[i-1] = current[i];
		i++;
	}
	delete [] current;
	return a;
}

/*Adds 'insert' row to end of vector*/
double* Vd_expandRow(double *current, int n1, int insert) {
	int i=0;
	double *a =Vd_alloc(n1+insert);
// Copy data before erase row
	for(i =0; i< n1; i++)
		a[i] = current[i];
	for(i =n1; i< n1+insert; i++)
		a[i] = 0;
	delete [] current;
	return a;
}

/* 2. DOUBLE MATRIX2D */
double **M2Dd_alloc(int n1, int n2) {
	double **a;
	a = new double* [sizeof(double*) * n1];
	a[0] = new double [sizeof(double) * n1 * n2];
	for (int i = 1; i < n1; i++)
		a[i] = a[0] + i * n2;
	return a;
}

double **M2Dd_resize(double **current, int row_old, int row_new, int col) {
	double **a;
	if(row_new == row_old)  return current;
	else if(row_new>row_old)
	{
		a=M2Dd_expandRow(current, row_old, col, row_new- row_old);
	}
	else
	{
		int i=0,j;
		a =M2Dd_alloc(row_new, col);
				// Copy data before erase row
		while(i<row_new)
		{
			for(j =0; j< col; j++)   a[i][j] = current[i][j];
			i++;
		}
//                std::cout << "Problem here2?" << std::endl;
		if(current!=NULL)  {delete [] current[0];   delete [] current;
		}
//                std::cout << "after2?" << std::endl;

	}
	return a;
}

/*Removes row form existing matrix*/
double ** M2Dd_eraseRow(double **current, int n1, int n2, int erase) {
	int i=0,j;
	double **a =M2Dd_alloc(n1-1, n2);
	// Copy data before erase row
	while(i<erase)
	{
		for(j =0; j< n2; j++)   a[i][j] = current[i][j];
		i++;
	}
	// Skip row = erase
	// Copy data after erased row
	i = erase+1;
	while(i<n1)
	{
		for(j =0; j< n2; j++)   a[i-1][j] = current[i][j];
		i++;
	}
	if(current!=NULL) {delete [] current[0];  delete [] current;
	}
	return a;
}

/*Adds 'insert' rows to end of matrix*/
double ** M2Dd_expandRow(double **current, int n1, int n2, int insert) {
	int i=0,j;
	double **a =M2Dd_alloc(n1+insert, n2);
// Copy data before erase row
	for(i =0; i< n1; i++)
		for(j =0; j< n2; j++)
		a[i][j] = current[i][j];
	for(i =n1; i< n1+insert; i++)
		for(j =0; j< n2; j++)
		a[i][j] = 0;

//   std::cout << "Problem here3?" << std::endl;
	if(current!=NULL)  {delete [] current[0];   delete [] current;}
//    std::cout << "after3?" << std::endl;

	return a;
}

//*******************************************************************************************/
// FLOAT MATRIX ALLOCATION:
//     5. FLOAT VECTOR
//     6. FLOAT MATRIX2D
//     7. FLOAT MATRIX3D
//     8. FLOAT MATRIX4D
//*******************************************************************************************/

/* 5. FLOAT VECTOR */
float *V_alloc(int n1) {
	float *V = new float[n1];
	return V;
}


/* 6. FLOAT MATRIX2D */
float **M2D_alloc(int n1, int n2) {
	float **a;
	a = new float* [sizeof(float*) * n1];
	a[0] = new float [sizeof(float) * n1 * n2];
	for (int i = 1; i < n1; i++)
		a[i] = a[0] + i * n2;
	return a;
}

/* -----------------------------------
INTEGER MATRIX ALLOCATION:
9. INTEGER VECTOR
10. INTEGER MATRIX2D
11. INTEGER MATRIX3D
12. INTEGER MATRIX4D
13. DOUBLE INTEGER VECTOR
----------------------------------- */
/* 9. INTEGER VECTOR */

int *Vi_alloc(int n1) {
	int *V = new int[n1];
	return V;
}


/* 10. INTEGER MATRIX2D */
int **IntM2D_alloc(int n1, int n2) {
	int **a;
	a = new int* [sizeof(int*) * n1];
	a[0] = new int [sizeof(int) * n1 * n2];
	for (int i = 1; i < n1; i++)
		a[i] = a[0] + i * n2;
	return a;
}


/* -----------------------------------
FREE MATRIX ALLOCATION:
14. FREE DOUBLE  VECTOR
15. FREE DOUBLE  MATRIX2D
16. FREE DOUBLE  MATRIX3D
17. FREE DOUBLE  MATRIX4D

18. FREE FLOAT   MATRIX2D
19. FREE FLOAT   MATRIX3D
20. FREE FLOAT   MATRIX4D

21. FREE INTEGER MATRIX2D
22. FREE INTEGER MATRIX3D
23. FREE INTEGER MATRIX4D
----------------------------------- */

/* 14. FREE VECTOR: use delete [] v  */

/* 15. FREE DOUBLE 2D. MATRIX */
void M2DdFree(double **matrix) {
	if(matrix!=NULL)
		{delete [] matrix[0];   delete [] matrix;}
}


/* 18. FREE FLOAT 2D. MATRIX */
void M2DFree(float **matrix) {
	if(matrix!=NULL)
	{
		delete [] matrix[0];   delete [] matrix;
		matrix=NULL;
	}
}

/* 21. FREE INTEGER 2D. MATRIX */
void M2DIFree(int **matrix) {
	if(matrix!=NULL)
	{
		delete [] matrix[0];   delete [] matrix;
		matrix=NULL;
	}
}


/*suggestion for Zap function

The delete and new operators in C++ are much better than the malloc and free functions of C. Consider using new and zap (delete function) instead of malloc and free as much as possible.

To make delete operators even more cleaner, make a Zap() inline function. Define a zap() function like this:

// Put an assert to check if x is NULL, this is to catch
// program "logic" errors early. Even though delete works
// fine with NULL by using assert you are actually catching
// "bad code" very early

// Defining Zap using templates
// Use zap instead of delete as this will be very clean
template <class T>
inline void zap(T & x) {
{assert(x != NULL);}
delete x;
x = NULL;
}

// In C++ the reason there are 2 forms of the delete operator is - because
// there is no way for C++ to tell the difference between a pointer to
// an object and a pointer to an array of objects. The delete operator
// relies on the programmer using "[]" to tell the two apart.
// Hence, we need to define zaparr function below.
// To delete array of pointers
template <class T>
inline void zaparr(T & x) {
{assert(x != NULL);}
delete [] x;
x = NULL;
}

The zap() function will delete the pointer and set it NULL. This will ensure that even if multiple zap()'s are called on the same deleted pointer then the program will not crash.
*/


// allocate a double vector_NR with subscript range v[nl..nh]
double *vector_NR(long nl, long nh) {
	double *v;
	v=(double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector_NR()");
	return v-nl+1;
}

// free a double vector_NR allocated with vector_NR()
void free_vector_NR(double *v, long nl, long nh) {
	free((char*) (v+nl-1));
}

// error handler
void nrerror(std::string error_text) {
	std::cerr << "run-time error...\n" << error_text<<std::endl;
}

