#include <cuda_runtime.h>

#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>
#include "TiMatrix.h"
#include "TiMatrixType.h"

using namespace std;

inline void CUDA_CALL(cudaError_t x)
{
	const cudaError_t a = x;
	if(a != cudaSuccess)
		cout<<"cuda error:"<<cudaGetErrorString(a)<<endl;

	//cudaDeviceReset();
	//assert(0);
}

__global__ void GetDiaKernel(CSR_Matrix DevA, float *D);
__global__ void DLMulKernel(CSR_Matrix DevK, CSR_Matrix DevA, float *D);
__global__ void CSRMulKernel(CSR_Matrix DevK, CSR_Matrix DevM, CSR_Matrix DevA);

float CSRMul(string fileName);


int main()
{
	float totaltime = 0;
	string filename;

	cout<<"please input matrix market file name:";
	cin>>filename;

	for(int i = 0 ; i < 5 ; i++)
	{
		totaltime += CSRMul(filename);
	}

	cout<<"the average cost time is "<<totaltime/5<<" ms "<<endl;

	return 0;

}



//·ÖÅädevice memory¸ø¾ØÕó
float CSRMul(string fileName)
{
	
	
	CSR_Matrix A;
	CSR_Matrix K;
	CSR_Matrix M;
	
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);
	
	ReadMatrixMarketFile(fileName.c_str(),A);
	cout<<"Reading matrix..."<<endl;

	//print the matrix
	cout<<A.n<<"    "<<A.n<<"    "<<A.nonzeroes<<endl;
	//cout<<"=================The Matrix A Format======================"<<endl;
	//PrintCSRFormat(A);

	
	
	K.n = A.n;
	K.nonzeroes = A.nonzeroes;
	K.mData = new float[A.nonzeroes];
	K.mIndex = new int[A.nonzeroes];
	K.mPtr = new int[A.n+1];


	M.n = A.n;
	M.nonzeroes = A.nonzeroes;
	M.mData = new float[A.nonzeroes];
	M.mIndex = new int[A.nonzeroes];
	M.mPtr = new int[A.n+1];

	CSR_Matrix DevK;
	CSR_Matrix DevM;
	CSR_Matrix DevA;
	float *DevD;

	float *D = new float[A.n];


	DevA.n = A.n;
	DevA.nonzeroes = A.nonzeroes;

	
	DevK.n = K.n;
	DevK.nonzeroes = K.nonzeroes;

	DevM.n = M.n;
	DevM.nonzeroes = M.nonzeroes;


	cudaMalloc((void **)&DevK.mData,sizeof(float)*A.nonzeroes);
	cudaMalloc((void **)&DevK.mIndex,sizeof(int)*A.nonzeroes);
	cudaMalloc((void **)&DevK.mPtr,sizeof(int)*(A.n+1));

	cudaMalloc((void **)&DevA.mData,sizeof(float)*A.nonzeroes);
	cudaMalloc((void **)&DevA.mIndex,sizeof(int)*A.nonzeroes);
	cudaMalloc((void **)&DevA.mPtr,sizeof(int)*(A.n+1));

	cudaMalloc((void **)&DevD, sizeof(float)*A.n);


	cudaMemcpy(DevA.mData,A.mData,sizeof(float)*A.nonzeroes,cudaMemcpyHostToDevice);
	cudaMemcpy(DevA.mIndex,A.mIndex,sizeof(int)*A.nonzeroes,cudaMemcpyHostToDevice);
	cudaMemcpy(DevA.mPtr,A.mPtr,sizeof(int)*(A.n+1),cudaMemcpyHostToDevice);

	cudaMemcpy(DevK.mPtr,A.mPtr,sizeof(int)*(A.n+1),cudaMemcpyHostToDevice);
	

	dim3 BlockSize = ceil(K.n/256.00);
	dim3 ThreadSize = 256;

	//Get the D from A
	//===================================

	cout<<"Entering the get D matrix kernel..."<<endl;
	GetDiaKernel<<<BlockSize,ThreadSize>>>(DevA, DevD);

	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
	CUDA_CALL(cudaMemcpy(D,DevD,sizeof(float)*(A.n),cudaMemcpyDeviceToHost));
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
	//cout<<"=================The Matrix D Format======================"<<endl;
	/*for(int i = 0 ; i < A.n ; i++)
		cout<<D[i]<<"    ";
	cout<<endl;*/
	


	//Get K from A
	cout<<"Entering the get K matrix kernel..."<<endl;
	//===================================
	DLMulKernel<<<BlockSize,ThreadSize>>>(DevK, DevA, DevD);

	printf("%s\n",cudaGetErrorString(cudaGetLastError()));

	cudaMemcpy(K.mData,DevK.mData,sizeof(float)*(A.nonzeroes),cudaMemcpyDeviceToHost);
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMemcpy(K.mIndex,DevK.mIndex,sizeof(int)*(A.nonzeroes),cudaMemcpyDeviceToHost);
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMemcpy(K.mPtr,DevK.mPtr,sizeof(int)*(A.n+1),cudaMemcpyDeviceToHost);
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));

	//cout<<"=================The Matrix K Format======================"<<endl;
	//PrintCSRFormat(K);


	cudaMalloc((void **)&DevM.mData,sizeof(float)*A.nonzeroes);
	cudaMalloc((void **)&DevM.mIndex,sizeof(int)*A.nonzeroes);
	cudaMalloc((void **)&DevM.mPtr,sizeof(int)*(A.n+1));

	cudaMemcpy(DevM.mPtr,A.mPtr,sizeof(int)*(A.n+1),cudaMemcpyHostToDevice);

	dim3 BlockSize1= ceil(A.n/8.0);
	dim3 ThreadSize1 = 256;


	//Get M from K
	//===================================
	cout<<"Entering the get M matrix kernel..."<<endl;
	CSRMulKernel<<<BlockSize1,ThreadSize1>>>(DevK,DevM,DevA);
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));

	cudaMemcpy(M.mData,DevM.mData,sizeof(float)*(M.nonzeroes),cudaMemcpyDeviceToHost);
	cudaMemcpy(M.mIndex,DevM.mIndex,sizeof(int)*(M.nonzeroes),cudaMemcpyDeviceToHost);
	cudaMemcpy(M.mPtr,DevM.mPtr,sizeof(int)*(M.n+1),cudaMemcpyDeviceToHost);

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	float elapseTime;
	cudaEventElapsedTime(&elapseTime,start,stop);

	cudaFree(DevK.mData);
	cudaFree(DevK.mIndex);
	cudaFree(DevK.mPtr);

	cudaFree(DevA.mData);
	cudaFree(DevA.mIndex);
	cudaFree(DevA.mPtr);

	cudaFree(DevM.mData);
	cudaFree(DevM.mIndex);
	cudaFree(DevM.mPtr);

	//cout<<"=================The Matrix M Format======================"<<endl;
	//PrintCSRFormat(M);
	delete [] A.mData;
	delete [] A.mIndex;
	delete [] A.mPtr;

	delete [] K.mData;
	delete [] K.mIndex;
	delete [] K.mPtr;

	delete [] M.mData;
	delete [] M.mIndex;
	delete [] M.mPtr;

	

	cout<<"the time cost is "<<elapseTime<<" ms."<<endl;
	return elapseTime;

}

__global__ void GetDiaKernel(CSR_Matrix DevA, float *D)
{
	int Row = blockDim.x * blockIdx.x+ threadIdx.x;
	int begin,end;

	if(Row < DevA.n)
	{

		begin = DevA.mPtr[Row];
		end = DevA.mPtr[Row+1];

		D[Row] = 0.0;

		
		for(int i = begin ; i < end ; i++)
		{
			if(Row == DevA.mIndex[i])
			{
				D[Row] = DevA.mData[i];
			}
		}

	}
	
}

__global__ void DLMulKernel(CSR_Matrix DevK, CSR_Matrix DevA, float *D)
{
	int Row = blockDim.x * blockIdx.x+ threadIdx.x;

	int begin,end;

	if(Row < DevA.n)
	{

		end = DevA.mPtr[Row+1];
		begin = DevA.mPtr[Row];

		for(int i = begin ; i < end ; i++)
		{
			int Col = DevA.mIndex[i];

			if(Row > Col)
			{
				DevK.mData[i] = -(DevA.mData[i]*(1/sqrt(D[Row]))*
					(1/D[Col]));
			}
			else if(Row == Col)
			{
				DevK.mData[i] = (1/sqrt(D[Row]));
			}
			else
				DevK.mData[i] = 0.0;

			DevK.mIndex[i] = DevA.mIndex[i];
		}
	}
}

__global__ void CSRMulKernel(CSR_Matrix DevK, CSR_Matrix DevM, CSR_Matrix DevA)
{
	//SpMM kernel for the CSR sparse matrix format using one 32-thread warp per matrix row of M

	__shared__ int RowIndex[8][1024]; 
	int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
	int warp_id = thread_id / 32;
	int lane = thread_id & (32 -1);

	int block_warp_id = warp_id %8;

	int Row = warp_id ;
	
	if(Row < DevA.n)
	{
		float dot = 0.0;
		int RowBegin = DevK.mPtr[Row];
		int RowEnd = DevK.mPtr[Row+1];

		for(int i = RowBegin ; i < RowEnd ; i++)
		{
			RowIndex[block_warp_id][i - RowBegin] = DevK.mIndex[i];
		}

		int beginA = DevA.mPtr[Row];
		int endA = DevA.mPtr[Row+1];

		for(int j = beginA + lane ; j < endA ; j += 32)
		{
			int Col = DevA.mIndex[j];

			int ColBegin = DevK.mPtr[Col];
			int ColEnd = DevK.mPtr[Col+1];

			for(int k = ColBegin ; k < ColEnd ; k++)
			{
				for(int l = RowBegin ; l < RowEnd ; l++)
				{
					if(DevK.mIndex[k] == RowIndex[block_warp_id][l-RowBegin])
					{
						dot += DevK.mData[k]*DevK.mData[l];
						break;
					}
				}
			}

			DevM.mData[j] += dot;
			DevM.mIndex[j] = DevA.mIndex[j];
		}
	}
}


