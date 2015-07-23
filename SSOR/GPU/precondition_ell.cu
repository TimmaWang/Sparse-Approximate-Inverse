#include <cuda_runtime.h>

#include <iostream>
#include <string>
#include <math.h>
#include "TiMatrix.h"
#include "TiMatrixType.h"

using namespace std;

__global__ void GetDiaKernel(ELL_Matrix DevA, float *D);
__global__ void DLMulKernel(ELL_Matrix DevK, ELL_Matrix DevA, float *D);
__global__ void ELLMulKernel(ELL_Matrix DevK, ELL_Matrix DevM, ELL_Matrix DevA);

float ELLMul(string filename);


int main()
{
	float totaltime = 0;
	string filename;

	cout<<"please input matrix market file name:";
	cin>>filename;

	for(int i = 0 ; i < 5 ; i++)
		totaltime += ELLMul(filename);

	cout<<"the average cost time is "<<totaltime/5<<" ms "<<endl;

	return 0;

}



//·ÖÅädevice memory¸ø¾ØÕó
float ELLMul(string filename)
{
	
	

	CSR_Matrix origin;
	ELL_Matrix A;
	ELL_Matrix K;
	ELL_Matrix M;
	float *D;
	
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);

	ReadMatrixMarketFile(filename.c_str(),origin);
	cout<<"Reading matrix..."<<endl;

	CSR2ELL(origin, A);

	

	//print the matrix
	cout<<A.nRow<<"    "<<A.nCol<<"    "<<A.nonzeroes<<endl;
	//PrintELLFormat(A);

	K.nRow = A.nRow;
	K.nCol = A.nCol;
	K.nonzeroes = A.nonzeroes;

	K.mData = new float[A.nRow*A.nCol];
	K.mIndex = new int[A.nRow*A.nCol];


	M.nRow = A.nRow;
	M.nonzeroes = A.nonzeroes;
	M.nCol = A.nCol;
	M.mData = new float[A.nRow*A.nCol];
	M.mIndex = new int[A.nRow*A.nCol];

	D = new float[A.nRow];

	ELL_Matrix DevK;
	ELL_Matrix DevM;
	ELL_Matrix DevA;
	float *DevD;


	DevA.nRow = A.nRow;
	DevA.nCol = A.nCol;
	DevA.nonzeroes = A.nonzeroes;

	DevK.nRow = K.nRow;
	DevK.nCol = K.nCol;
	DevK.nonzeroes = K.nonzeroes;

	DevM.nRow = M.nRow;
	DevM.nCol = M.nCol;
	DevM.nonzeroes = M.nonzeroes;


	cudaMalloc((void **)&DevK.mData,sizeof(float)*A.nRow*A.nCol);
	cudaMalloc((void **)&DevK.mIndex,sizeof(int)*A.nRow*A.nCol);

	cudaMalloc((void **)&DevA.mData,sizeof(float)*A.nRow*A.nCol);
	cudaMalloc((void **)&DevA.mIndex,sizeof(int)*A.nRow*A.nCol);

	cudaMalloc((void **)&DevD, sizeof(float)*A.nRow);


	cudaMemcpy(DevA.mData,A.mData,sizeof(float)*A.nRow*A.nCol,cudaMemcpyHostToDevice);
	cudaMemcpy(DevA.mIndex,A.mIndex,sizeof(int)*A.nRow*A.nCol,cudaMemcpyHostToDevice);
	

	dim3 BlockSize= ceil(K.nRow/256.00);
	dim3 ThreadSize = 256;

	cout<<"Entering the get D kernel..."<<endl;

	GetDiaKernel<<<BlockSize, ThreadSize>>>(DevA, DevD);
	cudaMemcpy(D,DevD,sizeof(float)*(A.nRow),cudaMemcpyDeviceToHost);

/*	for(int i = 0 ; i < A.nRow ; i++)
		cout<<D[i]<<" ";
	cout<<endl;*/

	//Get K from A
	cout<<"Entering the get K kernel..."<<endl;
	DLMulKernel<<<BlockSize,ThreadSize>>>(DevK, DevA, DevD);

	/*
	cudaMemcpy(K.mData,DevK.mData,sizeof(float)*(A.nRow*A.nCol),cudaMemcpyDeviceToHost);
	cudaMemcpy(K.mIndex,DevK.mIndex,sizeof(int)*(A.nRow*A.nCol),cudaMemcpyDeviceToHost);
	
	PrintELLFormat(K);
	*/
	cudaMalloc((void **)&DevM.mData,sizeof(double)*A.nRow*A.nCol);
	cudaMalloc((void **)&DevM.mIndex,sizeof(int)*A.nRow*A.nCol);



	dim3 BlockSize1= ceil(A.nRow/8.0);
	dim3 ThreadSize1 = 256;
	//Get M from K
	cout<<"Entering the get M kernel..."<<endl;
	ELLMulKernel<<<BlockSize1,ThreadSize1>>>(DevK,DevM,DevA);

	cudaMemcpy(M.mData,DevM.mData,sizeof(float)*(A.nRow*A.nCol),cudaMemcpyDeviceToHost);
	cudaMemcpy(M.mIndex,DevM.mIndex,sizeof(int)*(A.nRow*A.nCol),cudaMemcpyDeviceToHost);

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	float elapseTime;
	cudaEventElapsedTime(&elapseTime,start,stop);

	cout<<"the cost time is "<<elapseTime<<" ms "<<endl;

	cudaFree(DevK.mData);
	cudaFree(DevK.mIndex);

	cudaFree(DevA.mData);
	cudaFree(DevA.mIndex);

	cudaFree(DevM.mData);
	cudaFree(DevM.mIndex);

	cudaFree(DevD);

	//cout<<"=====================ELL M Matrix Format==========================="<<endl;
	//PrintELLFormat(M);

	delete [] A.mData;
	delete [] A.mIndex;

	delete [] K.mData;
	delete [] K.mIndex;

	delete [] M.mData;
	delete [] M.mIndex;

	delete [] D;

	return elapseTime;

}

__global__ void GetDiaKernel(ELL_Matrix DevA, float *D)
{
	int Row = blockDim.x * blockIdx.x+ threadIdx.x;
	int Col;

	if(Row < DevA.nRow)
	{
		D[Row] = 0.0;

		for(int i = 0 ; i < DevA.nCol ; i++)
		{
			Col = DevA.mIndex[i*DevA.nRow+Row];
			if(Row == Col)
			{
				D[Row] = DevA.mData[i*DevA.nRow+Row];
				break;
			}
		}
	}
}
__global__ void DLMulKernel(ELL_Matrix DevK, ELL_Matrix DevA, float *D)
{
	int Row = blockDim.x * blockIdx.x+ threadIdx.x;

	if(Row < DevA.nRow)
	{
		

		for(int i = 0 ; i < DevA.nCol ; i++)
		{
			int Col = DevA.mIndex[i*DevA.nRow+Row];

			if(Row > Col)
			{
				DevK.mData[i*DevA.nRow+Row] = -(DevA.mData[i*DevA.nRow+Row]*(1/sqrt(D[Row]))*
					(1/D[Col]));
			}
			else if(Row == Col)
			{
				DevK.mData[i*DevA.nRow+Row] = (1/sqrt(D[Row]));
			}
			else
				DevK.mData[i*DevA.nRow+Row] = 0.0;

			DevK.mIndex[i*DevA.nRow+Row] = DevA.mIndex[i*DevA.nRow+Row];
		}
		
	}
}
__global__ void ELLMulKernel(ELL_Matrix DevK, ELL_Matrix DevM, ELL_Matrix DevA)
{
	__shared__ int RowIndex[8][1024]; 
	int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
	int warp_id = thread_id / 32;
	int lane = thread_id & (32 -1);

	int block_warp_id = warp_id %8;

	int Row = warp_id ;

	if(Row < DevA.nRow)
	{
		float dot = 0.0;

		for(int i = 0 ; i < DevA.nCol ; i++)
		{
			RowIndex[block_warp_id][i] = DevK.mIndex[i*DevK.nRow+Row];
		}

		for(int l = lane ; l < DevA.nCol ; l+=32)
		{
			int Col = DevA.mIndex[l*DevK.nRow+Row];

			for(int i = 0 ; i < DevK.nCol ; i++)
			{
				int ColB = DevK.mIndex[i*DevK.nRow+Col];
			

				for(int j = 0 ; j < DevK.nCol ; j++)
				{
					int ColA = RowIndex[block_warp_id][j];
				
					if(ColA == ColB)
					{
						float valB =  DevK.mData[i*DevK.nRow+Col]; 
						float valA =  DevK.mData[j*DevK.nRow+Row];
							dot += valA*valB;
					}else if(ColB < ColA)
						break;
				}
			}
			DevM.mData[l*DevK.nRow+Row] += dot;
			DevM.mIndex[l*DevK.nRow+Row] = DevA.mIndex[l*DevK.nRow+Row];
		}
		
		//ELL*ELL
		
	}

}



