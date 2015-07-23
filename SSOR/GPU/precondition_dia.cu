#include <cuda_runtime.h>

#include <iostream>
#include <string>
#include <math.h>
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


__global__ void GetDiaKernel(DIA_Matrix DevA, float *D);
__global__ void DLMulKernel(DIA_Matrix DevK, DIA_Matrix DevA, float *D);
__global__ void DIAMulKernel(DIA_Matrix DevK, DIA_Matrix DevM, DIA_Matrix DevA);

float DIAMul(string filename);


int main()
{
	float totaltime = 0;
	string filename;

	cout<<"please input matrix market file name:";
	cin>>filename;

	for(int i = 0 ; i < 5 ; i++)
		totaltime += DIAMul(filename);

	cout<<"the average cost time is "<<totaltime/5<<" ms "<<endl;

	return 0;

}



//·ÖÅädevice memory¸ø¾ØÕó
float DIAMul(string filename)
{
	
	

	CSR_Matrix origin;
	DIA_Matrix A;
	DIA_Matrix K;
	DIA_Matrix M;
	float *D;
	
	ReadMatrixMarketFile(filename.c_str(),origin);
	cout<<"Reading matrix..."<<endl;

	CSR2DIA(origin, A);

	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);

	//print the matrix
	cout<<A.nRow<<"    "<<A.nCol<<"    "<<A.nonzeroes<<endl;
	//PrintDIAFormat(A);

	K.nRow = A.nRow;
	K.nCol = A.nCol;
	K.nonzeroes = A.nonzeroes;

	K.mData = new float[A.nRow*A.nCol];
	K.mOffset = new int[A.nCol];


	M.nRow = A.nRow;
	M.nCol = A.nCol;
	M.nDIG = A.nDIG;
	M.nonzeroes = A.nonzeroes;
	M.mData = new float[A.nRow*A.nCol];
	M.mOffset = new int[A.nCol];

	D = new float[A.nRow];

	DIA_Matrix DevK;
	DIA_Matrix DevM;
	DIA_Matrix DevA;
	float *DevD;


	DevA.nRow = A.nRow;
	DevA.nCol = A.nCol;
	DevA.nDIG = A.nDIG;
	DevA.nonzeroes = A.nonzeroes;

	DevK.nRow = K.nRow;
	DevK.nCol = K.nCol;
	DevK.nonzeroes = K.nonzeroes;

	DevM.nRow = M.nRow;
	DevM.nRow = M.nCol;
	DevM.nonzeroes = M.nonzeroes;


	cudaMalloc((void **)&DevK.mData,sizeof(float)*A.nRow*A.nCol);
	cudaMalloc((void **)&DevK.mOffset,sizeof(int)*A.nCol);

	cudaMalloc((void **)&DevA.mData,sizeof(float)*A.nRow*A.nCol);
	cudaMalloc((void **)&DevA.mOffset,sizeof(int)*A.nCol);

	cudaMalloc((void **)&DevD, sizeof(float)*A.nRow);


	cudaMemcpy(DevA.mData,A.mData,sizeof(float)*A.nRow*A.nCol,cudaMemcpyHostToDevice);
	cudaMemcpy(DevA.mOffset,A.mOffset,sizeof(int)*A.nCol,cudaMemcpyHostToDevice);

	cudaMemcpy(DevK.mOffset,A.mOffset,sizeof(int)*A.nCol,cudaMemcpyHostToDevice);
	

	dim3 BlockSize= ceil(K.nRow/256.00);
	dim3 ThreadSize = 256;

	cout<<"Entering the get D kernel..."<<endl;

	GetDiaKernel<<<BlockSize, ThreadSize>>>(DevA, DevD);
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMemcpy(D,DevD,sizeof(float)*(A.nRow),cudaMemcpyDeviceToHost);

	/*for(int i = 0 ; i < A.nRow ; i++)
		cout<<D[i]<<" ";
	cout<<endl;*/

	//Get K from A
	cout<<"Entering the get K kernel..."<<endl;
	DLMulKernel<<<BlockSize,ThreadSize>>>(DevK, DevA, DevD);
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
	
	//cout<<"=======================K Matirx===================="<<endl;
	cudaMemcpy(K.mData,DevK.mData,sizeof(float)*(A.nRow*A.nCol),cudaMemcpyDeviceToHost);
	cudaMemcpy(K.mOffset,DevK.mOffset,sizeof(int)*(A.nCol),cudaMemcpyDeviceToHost);
	//PrintDIAFormat(K);


	cudaMalloc((void **)&DevM.mData,sizeof(float)*A.nRow*A.nCol);
	cudaMalloc((void **)&DevM.mOffset,sizeof(int)*A.nCol);

	dim3 BlockSize1= ceil(A.nRow/8.0);
	dim3 ThreadSize1 = 256;
	//Get M from K
	cout<<"Entering the get M kernel..."<<endl;
	DIAMulKernel<<<BlockSize1,ThreadSize1>>>(DevK,DevM,DevA);
	printf("%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMemcpy(M.mData,DevM.mData,sizeof(float)*(A.nRow*A.nCol),cudaMemcpyDeviceToHost);
	cudaMemcpy(M.mOffset,DevM.mOffset,sizeof(int)*(A.nCol),cudaMemcpyDeviceToHost);

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	float elapseTime;
	cudaEventElapsedTime(&elapseTime,start,stop);

	cout<<"the cost time is "<<elapseTime<<" ms "<<endl;

	cudaFree(DevK.mData);
	cudaFree(DevK.mOffset);

	cudaFree(DevA.mData);
	cudaFree(DevA.mOffset);

	cudaFree(DevM.mData);
	cudaFree(DevM.mOffset);

	cudaFree(DevD);

	//cout<<"=======================M Matirx===================="<<endl;
	//PrintDIAFormat(M);

	delete [] A.mData;
	delete [] A.mOffset;

	delete [] K.mData;
	delete [] K.mOffset;

	delete [] M.mData;
	delete [] M.mOffset;


	return elapseTime;

}

__global__ void GetDiaKernel(DIA_Matrix DevA, float *D)
{
	int Row = blockDim.x * blockIdx.x+ threadIdx.x;
	int Col;

	if(Row < DevA.nRow)
	{
		D[Row] = 0.0;

		for(int i = 0 ; i < DevA.nCol ; i++)
		{
			Col = Row + DevA.mOffset[i];
			if(Row == Col)
			{
				D[Row] = DevA.mData[i*DevA.nRow+Row];
				break;
			}
		}
	}
}
__global__ void DLMulKernel(DIA_Matrix DevK, DIA_Matrix DevA, float *D)
{
	int Row = blockDim.x * blockIdx.x+ threadIdx.x;

	if(Row < DevA.nRow)
	{
		for(int i = 0 ; i < DevA.nCol ; i++)
		{
			int Col = Row + DevA.mOffset[i];

			if(Col >= 0 && Col < DevA.nRow)
			{
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

			}
		}
		
	}
}

__global__ void DIAMulKernel(DIA_Matrix DevK, DIA_Matrix DevM, DIA_Matrix DevA)
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
			RowIndex[block_warp_id][i] = Row + DevK.mOffset[i];
		}
		
		for(int l = lane ; l < DevA.nCol ; l+=32)
		{
			int Col = Row + DevA.mOffset[l];

			for(int i = 0 ; i < DevK.nCol ; i++)
			{
				int ColB = Col + DevK.mOffset[i];
			

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
			DevM.mOffset[l] = DevA.mOffset[l];
		}

		//DIA * DIA 
	}
}



