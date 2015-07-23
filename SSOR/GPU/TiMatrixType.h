#ifndef MTXFILE_H
#define MTXFILE_H

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;


typedef struct {
	int n ;	
	int nonzeroes ;
	int *mIndex ;
	int *mPtr ;
	float *mData ;
} CSR_Matrix ;

typedef struct {
	int nRow ;
	int nCol ;	
	int nDIG ;
	int nonzeroes ;
	int *mOffset ;
	float *mData ;
} DIA_Matrix ;

typedef struct {
	int nRow ;
	int nCol ;	
	int nonzeroes ;
	int *mIndex ;
	float *mData ;
} ELL_Matrix ;

typedef struct {
	float *vect ;
	int n   ;
} Vector ;

typedef struct {
	int n ;
	float *dgl ;
} DIAGS ;

#endif