#include "laplas_Solution.cuh"
#include <FL/Fl.H>
#include <fstream>
__global__ void setDouble1(double* target,double val){
	target[threadIdx.x+threadIdx.y*blockDim.x]=val;
}
__global__ void setDouble2(double* target,double val){
	target[threadIdx.x+threadIdx.y*blockDim.x+blockIdx.x*blockDim.x*blockDim.y]=/*(double)(threadIdx.x+threadIdx.y*blockDim.x+blockIdx.x*blockDim.x*blockDim.y)/(10*10*10)*/val;
}

///срез в кубической геометрии
__global__ void sliceKernel(double* A,double*result,const int z){
	int dY=blockDim.x;
	int dZ=(gridDim.x)*(blockDim.x);

	int i=threadIdx.x+blockIdx.x*dZ;
	double a=A[i+z*dY];

	result[threadIdx.x+blockIdx.x*dY]=a;
}

__global__ void laplasKernel(const double* fcu, double* fncu/*, double*sig,const double*Q*/){
	int dY=gridDim.x+2;
	int dZ=(gridDim.x+2)*(gridDim.y+2);
	double sig=0;
	int i0;
	i0=1+blockIdx.x +(blockIdx.y +1)*dY+(threadIdx.x +1)*dZ;
	int pts;
	//pts=4*M_PI;
	double upper_part;
	upper_part=0;
	double lower_part;
	lower_part=0;
	//int i,j,k;
	//for (i=-1;i<2;i+=2) {
	//	for(j=-1;j<2;j+=2){
	//		for(k=-1;k<2;k+=2){

	//			upper_part+=fcu[i0+i+j*dY+k*dZ]/sqrt(1.0*(i*i+j*j+k*k));
	//			lower_part+=1.0/sqrt(1.0*(i*i+j*j+k*k));

	//		}
	//	}
	//}
	//upper_part*=pts*sig;//[i0];
	//lower_part*=pts*sig;//[i0];
	upper_part+=1.0/6.0 * (fcu[i0+1]+fcu[i0-1]+fcu[i0+dY]+fcu[i0-dY]+fcu[i0+dZ]+fcu[i0-dZ]);
	lower_part+=1;
	fncu[i0]=upper_part/lower_part;

}
//__global__ void chargeKernel(
	__global__ void state_field_update(int* d_states,double* d_field_target){
		int dY=gridDim.x+2;
		int dZ=(gridDim.x+2)*(gridDim.y+2);
		int i0;
		i0=1+blockIdx.x +(blockIdx.y +1)*dY+(threadIdx.x +1)*dZ;
		if((d_states[i0]<30)&&(d_states[i0]>0)) d_field_target[i0]=1;
}

__global__ void yx_Borders(double* target){
	int dY=gridDim.x;
	int dZ=gridDim.x*blockDim.x;

	int idx=threadIdx.x+(blockIdx.x*dY);
	target[idx]=1;

	idx=threadIdx.x+(blockIdx.x*dY)+dZ*(dY-1);
	target[idx]=0;
}

__global__ void update_border(double* target){

	int dY=gridDim.x;
	int dZ=gridDim.x*blockDim.x;

	int idx=threadIdx.x+blockIdx.x*dZ; ///index of  y=const border
	target[idx]=target[idx+dY];/*
							   (double)(gridDim.x-blockIdx.x)/gridDim.x;*/
	__syncthreads();
	idx+=dZ-dY;
	target[idx]=target[idx-dY];/*
							   (double)(gridDim.x-blockIdx.x)/gridDim.x;*/
	__syncthreads();
	idx=threadIdx.x*dY+blockIdx.x*dZ;//x=const (0) border
	target[idx]=target[idx+1];/*
							  (double)(gridDim.x-blockIdx.x)/gridDim.x;*/
	__syncthreads();
	idx+=dY-1;// vtoraya x=const granica
	target[idx]=target[idx-1];/*
							  (double)(gridDim.x-blockIdx.x)/gridDim.x;*/
	__syncthreads();

}

__global__ void edge_update(double* target,int yx_idx,int dZ){
	int idx=yx_idx+threadIdx.x*dZ;   ///lesvie, stir'ek, igla ,elektrod
	target[idx]=1;
}

//// ABSOLUTE CONVERGENCE
__global__ void niv_stage_dim3(char* result,double* field_old,double* field_new, double eps){
	__syncthreads();
	int dY=gridDim.x;
	int dZ=gridDim.x*blockDim.x;
	int idx=threadIdx.x+blockIdx.x*dY+blockIdx.y*dZ;
	result[idx]= (abs(field_old[idx]-field_new[idx])<eps);
}
__device__ int niv_counter2;
__global__ void niv_stage_dim2(char* result){
	niv_counter2=0;
	__syncthreads();
	int dY=gridDim.x;
	int dZ=gridDim.x*blockDim.x;
	int idx=threadIdx.x+blockIdx.x*dY;
	for(int i=0;i<blockDim.x;i++)
		result[idx]=(result[idx]&&result[idx+i*dZ]);
	//atomicAdd(&niv_counter2,(int)result[idx]);
}
__device__ int niv_counter3;
__global__ void niv_stage_final(char*result){
	niv_counter3=0;
	int dY=gridDim.x;
	int idx=blockIdx.x;
	for(int i=0;i<blockDim.x;i++)
		result[idx]=(result[idx]&&result[idx+i*dY]);
	//atomicAdd(&niv_counter3, (int) result[idx]);
}

//// RELATIVE CONVERGENCE
__global__ void relNivStageDim3(char* result,double* field_old,double* field_new, double eps){
	__syncthreads();
	int dY=gridDim.x;
	int dZ=gridDim.x*blockDim.x;
	int idx=threadIdx.x+blockIdx.x*dY+blockIdx.y*dZ;
	result[idx]= (double)(abs(field_old[idx]-field_new[idx])/field_old[idx]<eps);
}

__global__ void relNivStageDim2(char* result){
	niv_counter2=0;
	__syncthreads();
	int dY=gridDim.x;
	int dZ=gridDim.x*blockDim.x;
	int idx=threadIdx.x+blockIdx.x*dY;
	for(int i=0;i<blockDim.x;i++)
		result[idx]=(result[idx]&&result[idx+i*dZ]);
	//atomicAdd(&niv_counter2,(int)result[idx]);
}

__global__ void relNivStageFinal(char*result){
	niv_counter3=0;
	int dY=gridDim.x;
	int idx=blockIdx.x;
	for(int i=0;i<blockDim.x;i++)
		result[idx]=(result[idx]&&result[idx+i*dY]);
	//atomicAdd(&niv_counter3, (int) result[idx]);
}

////
void cuLaplas::errReport(char* s, cudaError_t e){
	fprintf(stderr,s,e);
}

///инициализация массивов
cuLaplas::cuLaplas(int _size,int z,int limit,BORDER_CONDITION bc){
	///задаем чиселки
	_grSz=_size;
	_slc_z=z+1;
	dY=_grSz;
	dZ=_grSz*_grSz;
	///ограничили итерации
	_iteration_current=0;
	_iteration_total=0;
	_iteration_series=30;
	_iteration_limit=limit;
	///выбрали гран условие
	/*setBC(bc);*/
	errorsF= fopen("errlog.txt","w");
	///выбрали девайс
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(errorsF, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		//goto Error;
		return;
	}
	// Allocate GPU buffers()    .
	/////////////////////////////////////////////////////////////////////////////
	cudaStatus = cudaMalloc((void**)&dev_fi, _size *_size*_size* sizeof(double));
	cudaStatus = cudaMalloc((void**)&dev_niv_checks, _size *_size*_size* sizeof(char));
	cudaStatus = cudaMalloc((void**)&relDev_niv_checks, _size *_size*_size* sizeof(char));
	cudaStatus = cudaMalloc((void**)&dev_fi_old, _size *_size*_size* sizeof(double));
	cudaStatus = cudaMalloc((void**)&dev_fi_slice, _size *_size* sizeof(double));
	cudaStatus = cudaMalloc((void**)&devQ,_size*_size*_size*(sizeof(double)));
	cudaStatus = cudaMalloc((void**)&devQOld,_size*_size*_size*(sizeof(double)));
	////////////////////////////////////////////////////////////////////////////////

	////зануляем массивы
	////////////////////////////////////////////////////////////////////////////////

	cudaStatus=cudaMemset (dev_fi,0,sizeof(double)*_size*_size*_size);
	setDouble2
		<<<dim3(_size,1,1),dim3(_size,_size,1)>>>  (dev_fi,0.2);
	///////////////////////////////////
	cudaStatus=cudaMemset (dev_fi_old,0,sizeof(double)*_size*_size*_size);
	setDouble2
		<<<dim3(_size,1,1),dim3(_size,_size,1)>>>  (dev_fi_old,0.2);
	// Check for any errors
	//////////////////////////////////


	cudaStatus=cudaMemset (dev_fi_slice,0,sizeof(double)* _size * _size);
	cudaStatus=cudaMemset (dev_niv_checks,0,sizeof(char)* _size * _size*_size);
	cudaStatus=cudaMemset (relDev_niv_checks,0,sizeof(char)* _size * _size*_size);
	cudaStatus=cudaMemset (devQ,0,sizeof(float)*_size*_size*_size);
	cudaStatus=cudaMemset (devQOld,0,sizeof(float)*_size*_size*_size);
	//setDouble1<<<1,dim3(_size,_size,1)>>>(dev_fi_slice,0);
	//setDouble2<<<_size,dim3(_size,_size,1)>>>(dev_fi,0);
	//setDouble2<<<_size,dim3(_size,_size,1)>>>(dev_fi_old,0);
	cudaStatus=_neiman_init(dev_fi_old);
	cudaStatus=_neiman_init(dev_fi);
	//////////////////////////////////////////////////////////////////////////////////	
}

///срез хост
cudaError_t cuLaplas::cpySlice(double* host_target){
	sliceKernel<<<_grSz,_grSz>>>(dev_fi_old,dev_fi_slice,_slc_z);
	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.

	cudaStatus = cudaDeviceSynchronize();
	cudaStatus = cudaMemcpy(host_target, dev_fi_slice, _grSz *_grSz* sizeof(double), cudaMemcpyDeviceToHost);
	return cudaStatus;

}
cudaError_t cuLaplas::_neiman_noflow_Border(double* target,int current){

	update_border<<<_grSz,_grSz>>>(target);
	cudaStatus=cudaGetLastError();
	cudaStatus=cudaDeviceSynchronize();

	edge_update<<<1,_grSz/2>>>(target, _grSz*_grSz/2, _grSz*_grSz);
	cudaStatus=cudaGetLastError();
	cudaStatus=cudaDeviceSynchronize();
	return cudaStatus;
}
cudaError_t cuLaplas::_neiman_init(double* target){
	//////////
	yx_Borders<<<_grSz,_grSz>>>(target);
	cudaStatus=cudaGetLastError();
	cudaStatus=cudaDeviceSynchronize();
	/////////
	update_border<<<_grSz,_grSz>>>(target);
	cudaStatus=cudaGetLastError();
	cudaStatus=cudaDeviceSynchronize();
	///////////
	edge_update<<<1,_grSz/2>>>(target, _grSz*_grSz/2, _grSz*_grSz);
	cudaStatus=cudaGetLastError();
	cudaStatus=cudaDeviceSynchronize();
	/////////
	return cudaStatus;
}
cudaError_t cuLaplas::iteration(void* str_str,char* epsilon_check,double eps,float* time){

	float lapTime=0;
	cudaEvent_t start;
	cudaEvent_t stop;
	strmr_strct* a= (strmr_strct*)str_str;
	static int current=0;
	cudaEventCreate(&stop);
	cudaEventCreate(&start);
	cudaEventRecord(start,0);
	laplasKernel<<<dim3(_grSz-2,_grSz-2,1),dim3(_grSz-2,1,1)>>>(dev_fi_old,dev_fi/*,dev_sgm,dev_q*/);
	cudaStatus=cudaGetLastError();
	cudaStatus=cudaDeviceSynchronize();

	state_field_update<<<dim3(_grSz-2,_grSz-2,1),dim3(_grSz-2,1,1)>>>(a->_states,dev_fi/*,dev_sgm,dev_q*/);
	cudaStatus=cudaDeviceSynchronize();

	cudaStatus=_neiman_noflow_Border(dev_fi,current);

	//laplasKernel<<<dim3(_grSz-2,1,1),dim3(_grSz-2,_grSz-2,1)>>>(dev_fi,dev_fi_old/*,dev_sgm,dev_q*/);
	//cudaStatus=cudaGetLastError();
	//cudaStatus=cudaDeviceSynchronize();

	//cudaStatus=_neiman_noflow_Border(dev_fi_old,current);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&lapTime,start,stop);
	*time+=lapTime;
	current++;
	_iteration_total++;
	if (current>=_iteration_series){
		current=0;
		_iteration_current++;
	}
	swapFi();
	return cudaStatus;
}
void cuLaplas::convergence(char* epsilon_check, double eps){
	niv_stage_dim3<<<dim3(_grSz,_grSz,1),dim3(_grSz,1,1)>>>(dev_niv_checks,dev_fi_old,dev_fi,eps);
	cudaStatus=cudaDeviceSynchronize();

	niv_stage_dim2<<<dim3(_grSz,1,1),dim3(_grSz,1,1)>>>(dev_niv_checks);
	cudaStatus=cudaDeviceSynchronize();

	niv_stage_final<<<dim3(_grSz,1,1),1>>>(dev_niv_checks);
	cudaStatus=cudaDeviceSynchronize();

	cudaMemcpy(epsilon_check,dev_niv_checks,sizeof(char)*_grSz,cudaMemcpyDeviceToHost);
}
void cuLaplas::RelConvergence(char* epsilon_check, double eps){
	relNivStageDim3<<<dim3(_grSz,_grSz,1),dim3(_grSz,1,1)>>>(relDev_niv_checks,dev_fi_old,dev_fi,eps);
	cudaStatus=cudaDeviceSynchronize();

	relNivStageDim2<<<dim3(_grSz,1,1),dim3(_grSz,1,1)>>>(relDev_niv_checks);
	cudaStatus=cudaDeviceSynchronize();

	relNivStageFinal<<<dim3(_grSz,1,1),1>>>(relDev_niv_checks);
	cudaStatus=cudaDeviceSynchronize();

	cudaMemcpy(epsilon_check,relDev_niv_checks,sizeof(char)*_grSz,cudaMemcpyDeviceToHost);
}
cuLaplas:: ~cuLaplas(){
	cudaFree(dev_fi);
	cudaFree(dev_fi_old);
	cudaFree(dev_fi_slice);
}


///////////////////////////////////////////////////////////////////
///////////дальше идет код связанный с ростом стримерных каналов///
///////////////////////////////////////////////////////////////////

strmr_strct::strmr_strct(cuLaplas* parent){
	_size= parent->getSize();
	_slc_z=_size/2;
	cudaMalloc((void**)&_states,sizeof(int)*_size*_size*_size);
	cudaMemset (_states,0,sizeof(int)*_size*_size*_size);
	cudaMalloc((void**)&dev_slice,sizeof(int)*_size*_size);
	cudaMemset (dev_slice,0,sizeof(int)*_size*_size);
	cudaMalloc((void**)&rand_results,sizeof(float)*_size*_size*_size);
	initStates
		<<<dim3(_size-2,_size-2,1), dim3(_size-2,1,1)>>>
		(_states);

	edgeInitStruct
		<<<1,_size/2 -1>>>
		(_states, _size*_size/2, _size*_size);

	curandCreateGenerator(&gen, 
		CURAND_RNG_PSEUDO_DEFAULT);
}
//__device__ int gr_check(int* states,int id){
//	if( ) return 0;
//	return 1;
//}
cudaError_t strmr_strct::cpySlice(int* host_target){
	sliceKernel_int<<<_size,_size>>>(_states,dev_slice,_slc_z);
	// Check for any errors launching the kernel
	/*cudaStatus = cudaGetLastError();*/
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.

	/*cudaStatus = */cudaDeviceSynchronize();
	/*cudaStatus = */cudaMemcpy(host_target, dev_slice, _size *_size* sizeof(int), cudaMemcpyDeviceToHost);
	return cudaSuccess;

}
__global__ void initStates(int* A){
	int dY=gridDim.x+2;
	int dZ=(gridDim.x+2)*(gridDim.y+2);
	int i0=1+blockIdx.x +(blockIdx.y +1)*dY+(threadIdx.x +1)*dZ;
	A[i0]=200;

}
__global__ void edgeInitStruct(int*A,int xy_idx,int dZ){
	int i0=xy_idx+(1+threadIdx.x)*dZ;
	A[i0]=100;
	if ((threadIdx.x+1)==blockDim.x){
		A[i0]=101;
	}
}

///объявляем  счетчики
__device__ int d_count_200_true;
__device__ int d_count_101_true;
__device__ int d_count_0_30_true;
__device__ int d_count_31_60_true;
__device__ int d_count_calls;
__device__ int d_count_passed_checks[2];
__device__ void checkStructure(float rand_val,double* field,int* states,
							   int dY,int dZ,int id_to)
{
								   if	(states[id_to]>61){
									   if	(states[id_to]<=87){
										   /*atomicAdd(&d_count_31_60_true, 1);*/
										   int i=((states[id_to]-1)%9)%3-1;
										   int j=((states[id_to]-1)%9)/3-1;
										   int k=(states[id_to]-1)/9 -1;
										   int id_from =id_to-i-j*dY-k*dZ;
										   states[id_from]-=30;
										   states[id_to]-=30;
										   //atomicAdd(&d_count_31_60_true, 1);
									   }
								   }
								   if(states[id_to]==200 ){
									   /*atomicAdd(&d_count_200_true, 1);*/
									   for (int i=-1;i<2;i++){
										   for(int j=-1;j<2;j++){
											   for(int k=-1;k<2;k++){
												   int	id_from;
												   id_from=id_to+i+j*dY+k*dZ;
												   //if( states[id_from]<30){
												   //	if( states[id_from]>0){
												   //	if(!(i==0 && j==0 && k==0)){
												   //		__shared__ double randomized_field;
												   //		randomized_field=(field[id_from]-field[id_to])/sqrt((double)i*i+j*j+k*k)-log(rand_val)*0.1;
												   //		if (((i=1)|(i=-1))^
												   //			((j=1)|(j=-1))^
												   //			((k=1)|(k=-1)))
												   //		{   
												   //			if(1/*30.0/dY < randomized_field*/)						
												   //			states[id_to]=1+(i+1)+(j+1)*3+(k+1)*9;
												   //			return 1;
												   //		}
												   //		/*else {
												   //			if(30.0/dY < randomized_field)						
												   //			states[id_to]=31+(i+1)+(j+1)*3+(k+1)*9;
												   //			
												   //			return 1;
												   //		}*/
												   //	}	}
												   //	}

												   if( states[id_from]==101 || ((states[id_from]<60) &&(states[id_from]>0)))
												   {
													   /*atomicAdd(&d_count_101_true, 1);*/
													   if(!(i==0 && j==0 && k==0)){
														   double randomized_field;
														   randomized_field=abs(field[id_to]-field[id_from])/sqrt((double)i*i+j*j+k*k)-log(rand_val)*0.1;
														   if (((i==1)||(i==-1))^
															   ((j==1)||(j==-1))^
															   ((k==1)||(k==-1))){   
																   if(/*1000.0/dY*/0.90 < randomized_field){	
																	   if( (states[id_from]<60) &&(states[id_from]>30))
																	   {
																		   states[id_from]-=30;
																	   }
																	   states[id_to]=31+(i+1)+(j+1)*3+(k+1)*9;
																	   //atomicAdd(&d_count_passed_checks[0], 1);
																   }
														   }
														   else {
															   if(/*1.0/dY */0.90< randomized_field){						
																   states[id_to]=61+(i+1)+(j+1)*3+(k+1)*9;
																   //atomicAdd(&d_count_passed_checks[0], 1);
															   }
														   }

													   }
												   }
											   }
										   }
									   }
								   }
								   /*return 0;*/
}
__global__ void gr_iterate(int* states, double* field,float* uniformrand){
	d_count_200_true=0;
	d_count_101_true=0;
	d_count_0_30_true=0;
	d_count_31_60_true=0;

	int dY=gridDim.x+2;
	int dZ=(gridDim.x+2)*(gridDim.y+2);
	int id=1+blockIdx.x +(blockIdx.y +1)*dY+(threadIdx.x +1)*dZ;

	/*atomicAdd(&d_count_calls, 1);*/
	checkStructure(uniformrand[id],field,states,
		dY,dZ,id);

	__syncthreads();
}

__global__ void d_count_report(int* a){
	//printf("d_count_calls: %d  \n",d_count_calls);
	//printf("d_count_200_true: %d  \n",d_count_200_true);
	//printf("d_count_101_true: %d  \n",d_count_101_true);
	printf("d_count_passed_checks: %d  \n",d_count_passed_checks[0]);
	*a=d_count_passed_checks[0];
	/*
	printf("d_count_31_60_true: %d  \n",d_count_31_60_true);*/

	//printf("niv2: %d  \n",niv_counter2);
	//printf("niv3: %d  \n",niv_counter3);

}
int h_N_grow_counter=0;
void strmr_strct::count_report(/*std::fstream* a,*/double b){
	int* tmp;
	cudaMalloc((void**)&tmp,sizeof(int));
	d_count_report<<<1,1>>>(tmp);
	cudaDeviceSynchronize();
	cudaMemcpy(&h_N_grow_counter,tmp,sizeof(int),cudaMemcpyDeviceToHost);
	printf("h_N_grow_counter: %d  \n",h_N_grow_counter);
	/**a<<h_N_grow_counter<<'\n';*/
	printf("neodnorodnost: %f  \n",b);
}
__global__ void sliceKernel_int(int* A,int*result,const int z){
	int dY=blockDim.x;
	int dZ=(gridDim.x)*(blockDim.x);

	int idx=threadIdx.x+blockIdx.x*dZ;

	int a;
	for(int i=0;i<gridDim.x;i++){
		a=A[idx+i*dY];
		if((a!=200)&&(a!=0))break;
	}
	result[threadIdx.x+blockIdx.x*dY]=a;
}

void strmr_strct::cu_iterate(double* field){
	curandSetPseudoRandomGeneratorSeed(gen,seed_rng());
	curandGenerateUniform(gen,rand_results,(_size)*(_size)*(_size));
	cudaDeviceSynchronize();
	gr_iterate<<<dim3(_size-2,_size-2,1),
		dim3(_size-2,1,1)>>>
		(_states, field, rand_results);
	cudaDeviceSynchronize();


}
