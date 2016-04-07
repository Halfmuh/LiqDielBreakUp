
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "curand.h"
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <cstdio>

//#include <cstdint>
#include <stdio.h>

//__constant__ int *dev_size;
//__constant__ int *dev_slice_z;
enum BORDER_CONDITION{
	NOFLOW,
	LINEAR
};
///////////////////////////////////////////////////////////////
///////для расчета электрических полей/////////////////////////
///////////////////////////////////////////////////////////////
class cuLaplas{
	int _grSz;
	int dY;
	int dZ;
	int _slc_z;
	int _iteration_current; ///текущая серия
	int _iteration_series;///размер серии итераций
	int _iteration_total;///общее количество завершенных итераций
	int _iteration_limit;  ////предельное количество серий итераций
	//cudaError_t (*cuLaplas::current_bc) (double*);
	cudaError_t _neiman_noflow_Border(double*,int current);
	cudaError_t _neiman_init(double*);
	/*void setBC(BORDER_CONDITION bc){
	switch (bc){
	case NOFLOW:
	this->current_bc= _neiman_noflow_Border;
	break;
	}
	}*/
public:
	FILE* errorsF;
	cuLaplas(int _size,int z,int limit,BORDER_CONDITION bc);////
	double* dev_fi;
	double* dev_fi_old;
	double* dev_fi_slice;

	char* dev_niv_checks;
	char* relDev_niv_checks;
	//double* dev_sgm;
	//float* devQ;
	//float* devQOld;
	cudaError_t cudaStatus;
	int getTotal(){
		return _iteration_total;
	}
	int getCurrent(){
		return _iteration_current;
	}
	int getSize(){
		return _grSz;
	}
	void setSlice(int _slc){
		_slc_z=_slc;
	}
	void errReport(char* s, cudaError_t);
	cudaError_t cpySlice(double* target);
	cudaError_t iteration(void* str_str,char* epsilon_check,double eps,float* time);
	void convergence(char* epsilon_check, double eps);
	void RelConvergence(char* epsilon_check, double eps);
	void swapFi(){
		double* tmp= dev_fi;
		dev_fi=dev_fi_old;
		dev_fi_old=tmp;
	}
	~cuLaplas();
};
///////////////////////////////////////////////////////////////
///////для роста каналов///////////////////////////////////////
///////////////////////////////////////////////////////////////
class strmr_strct{
	int _size;
public:
	std::random_device seed_rng;
	curandGenerator_t gen;
	int* _states;
	float* rand_results;
	int* dev_slice;
	int _slc_z;
	///////////////////////////
	///0-граница
	///1-27 - стримеры
	///31-57 -будут стримерами
	/// 100 -острие
	///200 - диэлектрик
	///////////////////////////
	strmr_strct( cuLaplas* parent);
	void cu_iterate(double* field);
	void count_report(/*std::fstream* a,*/double b);
	cudaError_t cpySlice(int* host_target);
	~strmr_strct();
};


///////////////////////////////////////////////////////////////
/////////////расчет полей на гпу//////////////////////////////
//////////////////////////////////////////////////////////////
__global__ void sliceKernel(double* A,
							double*result,
							const int x);
__global__ void laplasKernel(const double* fcu,
							 double* fncu/*,double*sig,const double*Q*/);
__global__ void yx_Borders(double* target);
__global__ void update_border(double* target);
__global__ void edge_update(double* target,
							int yx_idx,
							int dZ);


///////////////////////////////////////////////////////////////
////////////рост структуры на гпу//////////////////////////////
///////////////////////////////////////////////////////////////
__global__ void sliceKernel_int(int* A,int*result,const int z);
__global__ void initStates(int* A);
__global__ void edgeInitStruct(int*A,int xy_idx,int dZ);
__device__ int checkStructure(float rand_val,double* field,int states,
							  int dY,int dZ,int id_to);
__device__ int strmr_update();
__global__ void gr_iterate(int* satates, double* field,float* uniformrand);

__global__ void strmr_growth();
