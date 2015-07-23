#ifndef MATRIXOPERATE_H
#define MATRIXOPERATE_H

#include "TiMatrixType.h"
#include <iomanip>
#include <stdlib.h>


//从高到低排序，快排的参数
int compInc(const void *a, const void *b)  
{  
  return *(int*)a - *(int*)b;  
}

//判断下标是否存在
bool isExist(int *diaIndex, int diaNum, int newIndex)
{
	bool rValue = false ;
	for( int i = 0 ; i < diaNum ; i++ ){
		if( diaIndex[i] == newIndex ){
			rValue = true ;
			break ;
		}
	}
	return rValue ;
}

//寻找下标的位置
int findIndexPos(int *diaIndex, int nCol, int newIndex)
{
	int rValue = nCol + 2 ;
	for( int i = 0 ; i < nCol ; i++ ){
		if( diaIndex[i] == newIndex ){
			rValue = i ;
			break ;
		}
	}
	return rValue ;
}

//从矩阵市场获得的矩阵文件.mtx，将其中的矩阵按照CSC的形式存储
void ReadMatrixMarketFileCSC(const char *fileName,CSR_Matrix &A)
{
	ifstream read;

	int currentRow = 0;
	int nextRow = 0;
	int num_per_row = 0;
	int count = 0;

	read.open(fileName);

	if(!read)
	{
		cout<<"can not open the file"<<endl;
		exit(1);
	}

	read>>A.n;
	read>>A.n;
	read>>A.nonzeroes;
	
	A.mData = new float[A.nonzeroes];
	A.mIndex = new int[A.nonzeroes];
	A.mPtr = new int[A.n+1];

	A.mPtr[0] = 0;
	int i ;
	int j ;
	float data ;

	while(read >> i
		&&read >> j
		&&read >> data)
	{
		
		i = i-1;
		j = j-1;
		//cout<<i<<" "<<j<<" "<<data<<endl;
		if(currentRow == j)
		{
			//num_per_row++;
			A.mData[count] = data;
			A.mIndex[count] = i;

		}
		else
		{
			nextRow = j;
			int rowInter = nextRow - currentRow;

			A.mPtr[currentRow+1] = count;

			while(--rowInter)
			{
				A.mPtr[++currentRow] = count;
			}

			currentRow = j;
			A.mData[count] = data;
			A.mIndex[count] = i;
		}
		
		//read.getline(temp,100);
		count++;
	}

	A.mPtr[currentRow+1] = count; 
}
//从矩阵市场获得的矩阵文件.mtx，将其中的矩阵按照CSR的形式存储
void ReadMatrixMarketFile(const char *fileName,CSR_Matrix &A)
{
	/*
	ifstream read;

	int currentRow = 0;
	int nextRow = 0;
	int num_per_row = 0;
	int count = 0;

	read.open(fileName);

	if(!read)
	{
		cout<<"can not open the file"<<endl;
		exit(1);
	}

	char temp[100];
	read.getline(temp,100);
	
	string s = temp;
	string s1 = s.substr(0,s.find_first_of(' '));

	string s2 = s.substr(s.find_first_of(' ')+1);
	string s3 = s2.substr(0,s.find_first_of(' '));

	string s4 = s2.substr(s2.find_first_of(' ')+1);

	A.n = atoi(s1.c_str());
	A.n = atoi(s3.c_str());
	A.nonzeroes = atoi(s4.c_str());

	A.mData = new float[A.nonzeroes];
	A.mIndex = new int[A.nonzeroes];
	A.mPtr = new int[A.n+1];

	A.mPtr[0] = 0;

	
	read.getline(temp,100);
	while(!read.eof())
	{
		s = temp;
		s1 = s.substr(0,s.find_first_of(' '));

		s2 = s.substr(s.find_first_of(' ')+1);
		s3 = s2.substr(0,s.find_first_of(' '));

		s4 = s2.substr(s2.find_first_of(' ')+1);

		int i = atoi(s1.c_str()) -1;
		int j = atoi(s3.c_str()) -1;
		float data = atof(s4.c_str());

		if(currentRow == i)
		{
			//num_per_row++;
			A.mData[count] = data;
			A.mIndex[count] = j;

		}
		else
		{
			nextRow = i;
			int rowInter = nextRow - currentRow;

			A.mPtr[currentRow+1] = count;

			while(--rowInter)
			{
				A.mPtr[++currentRow] = count;
			}

			currentRow = i;
			A.mData[count] = data;
			A.mIndex[count] = j;
		}
		
		read.getline(temp,100);
		count++;
	}

	A.mPtr[currentRow+1] = count; 
	*/

	ifstream read;

	int currentRow = 0;
	int nextRow = 0;
	int num_per_row = 0;
	int count = 0;

	read.open(fileName);

	if(!read)
	{
		cout<<"can not open the file"<<endl;
		exit(1);
	}

	read>>A.n;
	read>>A.n;
	read>>A.nonzeroes;
	
	A.mData = new float[A.nonzeroes];
	A.mIndex = new int[A.nonzeroes];
	A.mPtr = new int[A.n+1];

	A.mPtr[0] = 0;
	int i ;
	int j ;
	float data ;

	while(read >> i
		&&read >> j
		&&read >> data)
	{
		
		i = i-1;
		j = j-1;
		//cout<<i<<" "<<j<<" "<<data<<endl;
		if(currentRow == i)
		{
			//num_per_row++;
			A.mData[count] = data;
			A.mIndex[count] = j;

		}
		else
		{
			nextRow = i;
			int rowInter = nextRow - currentRow;

			A.mPtr[currentRow+1] = count;

			while(--rowInter)
			{
				A.mPtr[++currentRow] = count;
			}

			currentRow = i;
			A.mData[count] = data;
			A.mIndex[count] = j;
		}
		
		//read.getline(temp,100);
		count++;
	}

	A.mPtr[currentRow+1] = count; 
}

//Get ELL Format From CSR Format
void CSR2ELL(const CSR_Matrix CSR_A, ELL_Matrix &ELL_A)
{
	int nRow = CSR_A.n ;
	
	int maxBand = 0 ;
	for( int i = 0 ; i < nRow ; i++ ){
		int nRowStart = CSR_A.mPtr[i] ;
		int nRowEnd = CSR_A.mPtr[i+1] ;
		if( maxBand < (nRowEnd - nRowStart) ){
			maxBand = nRowEnd - nRowStart ;
		}
	}
	
	ELL_A.nRow = nRow ;
	ELL_A.nCol = maxBand ;
	ELL_A.nonzeroes = CSR_A.nonzeroes ;
	
	ELL_A.mData = ( float* )calloc(  ELL_A.nRow * ELL_A.nCol, sizeof( float ) ) ;
	ELL_A.mIndex = ( int* )calloc(  ELL_A.nRow * ELL_A.nCol, sizeof( int ) ) ;
	
	//compute the ELL_A->mData and ELL_A->mIndex
	for( int i = 0 ; i < nRow ; i++ ){
		int nRowStart = CSR_A.mPtr[i] ;
		int nRowEnd = CSR_A.mPtr[i+1] ;
		for( int j = nRowStart ; j < nRowEnd ; j++ ){
			ELL_A.mData[i + (j - nRowStart)*nRow] = CSR_A.mData[j] ;
			ELL_A.mIndex[i + (j - nRowStart)*nRow] = CSR_A.mIndex[j] ;
		}
	}
	
}


//Get DIA Format from CSR Format
void CSR2DIA(const CSR_Matrix CSR_A, DIA_Matrix &DIA_A)
{
	
	int nRow = CSR_A.n ;
	
	int diaNum = 0 ;
	int *diaIndex ; 
  
	diaIndex = (int*)malloc( sizeof(int) * (nRow * 2 + 1) ) ;
  
  //--compute the number of diagonal lines and offsets
	for( int i = 0 ; i < nRow ; i++ ){
		int nRowStart = CSR_A.mPtr[i] ;
		int nRowEnd = CSR_A.mPtr[i+1] ;
		for( int j = nRowStart ; j < nRowEnd ; j++ ){
			if(diaNum == 0){
				diaIndex[diaNum] = CSR_A.mIndex[j] - i ;
				diaNum++ ;
			}else{
			  int newIndex = CSR_A.mIndex[j] - i ;
			  if( !isExist(diaIndex, diaNum, newIndex) ){
			  	diaIndex[diaNum] = newIndex ;
			  	diaNum++ ;
			  }
			}
		}
	}
	
	//sort elements of diaIndex from low to high
	qsort(diaIndex, diaNum, sizeof(diaIndex[0]), compInc);
	
	//--copy diaIndex and diaNum to DIA_A->mOffset and DIA_A->nCol
	DIA_A.nCol = diaNum ;
	DIA_A.mOffset = (int*)malloc( sizeof( int ) * diaNum ) ;
	
	for( int i = 0 ; i < diaNum ; i++ ){
		DIA_A.mOffset[i] = diaIndex[i] ;
	}
	
	free( diaIndex ) ;
	
	DIA_A.nDIG = nRow ;
	DIA_A.nRow = nRow ;
	DIA_A.nonzeroes = CSR_A.nonzeroes ;
	
	DIA_A.mData = (float*)calloc( DIA_A.nRow * DIA_A.nCol, sizeof( float )  ) ;
	
	//--compute the DIA_A->mData 
	for( int i = 0 ; i < nRow ; i++ ){
		int nRowStart = CSR_A.mPtr[i] ;
		int nRowEnd = CSR_A.mPtr[i+1] ;
		for( int j = nRowStart ; j < nRowEnd ; j++ ){
			int newIndex = CSR_A.mIndex[j] - i ;
			int iPos = findIndexPos( DIA_A.mOffset, DIA_A.nCol, newIndex ) ;
			DIA_A.mData[i + iPos * nRow ] = CSR_A.mData[j] ;
		}
	} 
}

void PrintCSRFormat(const CSR_Matrix A)
{
	cout<<A.n<<"    "<<A.n<<"    "<<A.nonzeroes<<endl;
	for(int i = 0 ; i < A.nonzeroes ; i++)
		cout<<setprecision(5)<<A.mData[i]<<"   ";
	cout<<endl;
	cout<<"================================================="<<endl;
	for(int i = 0 ; i < A.nonzeroes ; i++)
		cout<<A.mIndex[i]<<"   ";
	cout<<endl;
	cout<<"================================================="<<endl;

	for(int i = 0 ; i < A.n+1 ; i++)
		cout<<A.mPtr[i]<<"   ";
	cout<<endl;
	cout<<"================================================="<<endl;
}



void PrintDIAFormat(const DIA_Matrix A)
{
	cout<<A.nRow<<"    "<<A.nCol<<"    "<<A.nonzeroes<<endl;
	cout<<A.nDIG<<endl;
	cout<<"==============A Data=========================="<<endl;
	for(int i = 0 ; i < A.nRow ; i++)
	{
		for(int j = i ; j < A.nCol*A.nRow ; j+=A.nRow)
			cout<<A.mData[j]<<" ";
		cout<<endl;
	}
	cout<<endl;
	cout<<"=================================================="<<endl;
	cout<<"==============A Offset=========================="<<endl;
	for(int i = 0 ; i < A.nCol ; i++)
	{
		cout<<setw(4)<<left<<setprecision(3)<<A.mOffset[i];
	}
	cout<<endl;
	cout<<"=================================================="<<endl;
}


void PrintELLFormat(const ELL_Matrix A)
{
	cout<<A.nRow<<"    "<<A.nCol<<"    "<<A.nonzeroes<<endl;
	cout<<"===============A Data========================="<<endl;
	for(int i = 0 ; i < A.nRow ; i++)
	{
		for(int j = i ; j < A.nCol*A.nRow ; j+=A.nRow)
			cout<<A.mData[j]<<" ";
		cout<<endl;
	}
	cout<<endl;
	cout<<"========================================"<<endl;
	cout<<"===============A Index========================="<<endl;
	for(int i = 0 ; i < A.nRow ; i++)
	{
		for(int j = i ; j < A.nCol*A.nRow ; j+=A.nRow)
			cout<<A.mIndex[j]<<" ";
		cout<<endl;
	}
	cout<<endl;
	cout<<"========================================"<<endl;
}

/******************************************************
function name : CSR2CSRPadding
parameter     : CSR_Matrix origin, CSR_Matrix padding
return        : void
*******************************************************/
void CSR2CSRPadding(const CSR_Matrix origin, CSR_Matrix &padding)
{
	int begin,end,row_num,sum;
	int beginP,endP;
	sum =0;
	padding.mPtr = (int*)calloc( origin.n+1, sizeof( int )  ) ; 
	padding.mPtr[0] = 0;
	//compute each row number of origin matrix
	for(int i = 0 ; i < origin.n ; i++)
	{
		begin = origin.mPtr[i];
		end = origin.mPtr[i+1];

		row_num = ((end - begin)/32 + 1)*32;
		sum += row_num;

		padding.mPtr[i+1] = sum;

	}

	padding.n = origin.n;
	padding.nonzeroes = sum;
	padding.mData = (float*)calloc( sum, sizeof( float )  ) ; 
	padding.mIndex = (int*)calloc( sum, sizeof( int )  ) ; 

	//copy the origin elements to padding
	for(int i = 0 ; i < padding.n ; i++)
	{
		int count = 0 ;
		beginP = padding.mPtr[i];
		endP = padding.mPtr[i+1];

		begin = origin.mPtr[i];
		end = origin.mPtr[i+1];

		while(begin+count < end)
		{
			padding.mData[beginP+count] = origin.mData[begin+count];
			padding.mIndex[beginP+count] = origin.mIndex[begin+count];
			count++;
		}

	}
}

/******************************************************
function name : CSR2ELLPadding
parameter     : const CSR_Matrix CSR_A, ELL_Matrix &ELL_A
return        : void
*******************************************************/
void CSR2ELLPadding(const CSR_Matrix CSR_A, ELL_Matrix &ELL_A)
{
	int nRow = CSR_A.n ;
	
	int maxBand = 0 ;
	for( int i = 0 ; i < nRow ; i++ ){
		int nRowStart = CSR_A.mPtr[i] ;
		int nRowEnd = CSR_A.mPtr[i+1] ;
		if( maxBand < (nRowEnd - nRowStart) ){
			maxBand = nRowEnd - nRowStart ;
		}
	}
	
	ELL_A.nRow = nRow ;
	//padding the nCol
	ELL_A.nCol = (maxBand /32 +1)*32;
	ELL_A.nonzeroes = CSR_A.nonzeroes ;
	
	ELL_A.mData = ( float* )calloc(  ELL_A.nRow * ELL_A.nCol, sizeof( float ) ) ;
	ELL_A.mIndex = ( int* )calloc(  ELL_A.nRow * ELL_A.nCol, sizeof( int ) ) ;
	
	//compute the ELL_A->mData and ELL_A->mIndex
	for( int i = 0 ; i < nRow ; i++ ){
		int nRowStart = CSR_A.mPtr[i] ;
		int nRowEnd = CSR_A.mPtr[i+1] ;
		for( int j = nRowStart ; j < nRowEnd ; j++ ){
			ELL_A.mData[i + (j - nRowStart)*nRow] = CSR_A.mData[j] ;
			ELL_A.mIndex[i + (j - nRowStart)*nRow] = CSR_A.mIndex[j] ;
		}
	}
}


/******************************************************
function name : CSR2DIAPadding
parameter     : const CSR_Matrix CSR_A, DIA_Matrix &DIA_A
return        : void
*******************************************************/
void CSR2DIAPadding(const CSR_Matrix CSR_A, DIA_Matrix &DIA_A)
{
	int nRow = CSR_A.n ;
	
	int diaNum = 0 ;
	int *diaIndex ; 
  
	diaIndex = (int*)malloc( sizeof(int) * (nRow * 2 + 1) ) ;
  
  //--compute the number of diagonal lines and offsets
	for( int i = 0 ; i < nRow ; i++ ){
		int nRowStart = CSR_A.mPtr[i] ;
		int nRowEnd = CSR_A.mPtr[i+1] ;
		for( int j = nRowStart ; j < nRowEnd ; j++ ){
			if(diaNum == 0){
				diaIndex[diaNum] = CSR_A.mIndex[j] - i ;
				diaNum++ ;
			}else{
			  int newIndex = CSR_A.mIndex[j] - i ;
			  if( !isExist(diaIndex, diaNum, newIndex) ){
			  	diaIndex[diaNum] = newIndex ;
			  	diaNum++ ;
			  }
			}
		}
	}
	
	//sort elements of diaIndex from low to high
	qsort(diaIndex, diaNum, sizeof(diaIndex[0]), compInc);
	
	//--copy diaIndex and diaNum to DIA_A->mOffset and DIA_A->nCol
	//padding the nCol
	DIA_A.nCol = (diaNum/32 + 1)*32 ;
	DIA_A.mOffset = (int*)malloc( sizeof( int ) * diaNum ) ;
	
	for( int i = 0 ; i < diaNum ; i++ ){
		DIA_A.mOffset[i] = diaIndex[i] ;
	}
	
	free( diaIndex ) ;
	
	DIA_A.nDIG = nRow ;
	DIA_A.nRow = nRow ;
	DIA_A.nonzeroes = CSR_A.nonzeroes ;
	
	DIA_A.mData = (float*)calloc( DIA_A.nRow * DIA_A.nCol, sizeof( float )  ) ;
	
	//--compute the DIA_A->mData 
	for( int i = 0 ; i < nRow ; i++ ){
		int nRowStart = CSR_A.mPtr[i] ;
		int nRowEnd = CSR_A.mPtr[i+1] ;
		for( int j = nRowStart ; j < nRowEnd ; j++ ){
			int newIndex = CSR_A.mIndex[j] - i ;
			int iPos = findIndexPos( DIA_A.mOffset, DIA_A.nCol, newIndex ) ;
			DIA_A.mData[i + iPos * nRow ] = CSR_A.mData[j] ;
		}
	} 
}

#endif