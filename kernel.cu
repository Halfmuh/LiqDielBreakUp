#if HAVE_PTHREAD || defined(WIN32)
#include "threads.h"
#include "Control_window.h"
#include "Plot_Window.h"
#include "laplas_Solution.cuh"
#include <iostream>
#include <cuda_profiler_api.h>
#include <fstream>
#include <sstream>
//#include <cstdint>
////оконца
Control_Window *ctrlwnd;
Plot_Window *pltwnd;
Plot_Window *pltwnd2;
////поток
Fl_Thread fl_thread_draw;
Fl_Thread fl_thread_calc;
///массивы

cuLaplas* lapls;
strmr_strct* str_str;
int* transfer_storage;
double* output;



static void start_cb(Fl_Widget *w,void* data);
static void stop_cb(Fl_Widget *w,void* data);
void thread_calc(void* p);///это основна€ функци€ работ€юща€ в потоке расчета
void thread_draw(void* p);

///// save to file binary data
///// length - array size

bool saveArray( const double* pdata, size_t length, const std::string& file_path )
{
	std::ofstream os(file_path.c_str(), std::ios::binary | std::ios::out);
	if ( !os.is_open() )
		return false;
	os.write(reinterpret_cast<const char*>(pdata), std::streamsize(length*sizeof(double)));
	os.close();
	return true;
}

bool loadArray( double* pdata, size_t length, const std::string& file_path)
{
	std::ifstream is(file_path.c_str(), std::ios::binary | std::ios::in);
	if ( !is.is_open() )
		return false;
	is.read(reinterpret_cast<char*>(pdata), std::streamsize(length*sizeof(double)));
	is.close();
	return true;
}

/////
/////

int main()
{
	int count = 0;
	cudaGetDeviceCount(&count);
	ctrlwnd=new Control_Window(200,400,start_cb,stop_cb,"ctrl") ;

	((Fl_Window*)ctrlwnd)->show();
	while(Fl::wait()>0){
		if((Plot_Window*)Fl::thread_message()==pltwnd){

			if(pltwnd!=0){
				Fl::lock();
				printf("redraw called\n");
				//for (int i=0;i<pltwnd->getGridSize()*pltwnd->getGridSize();i++)pltwnd->arr[i]=transfer_storage[i];
				//for(int i=0;i<pltwnd->getGridSize();i++) for(int j=0;j<pltwnd->getGridSize();j++)
				//	std::cout<<transfer_storage[i+j*pltwnd->getGridSize()]<<" ";
				for(int i=0;i<pltwnd2->getGridSize();i++){
					for(int j=0;j<pltwnd2->getGridSize();j++){
						pltwnd2->arr[i+j*pltwnd->getGridSize()]=transfer_storage[i+j*pltwnd->getGridSize()];
					}}
				pltwnd->redraw();
				pltwnd2->redraw();
				Fl::unlock();
			}
		}
		else if((cuLaplas*)Fl::thread_message()==lapls)
			if(lapls!=0) lapls->errReport("error code %d",lapls->cudaStatus);
	}
	return cudaDeviceReset();
}



///////////To Do  разгрести эти callback'и
void start_cb(Fl_Widget *w,void* data){
	Control_Window* tmp =(Control_Window*)w->parent();

	if(tmp->aftertasks & 2)
	{
		delete pltwnd;
		tmp->aftertasks^=2;
	}

	if((tmp->calc_state) == 0)
	{
		tmp->calc_state =1;
		///создадим оконце
		pltwnd= new Plot_Window(900,900,tmp,"Plot");
		pltwnd2= new Plot_Window(900,900,tmp,"Plot2");
		///и заведм массивы
		lapls=new cuLaplas(pltwnd->getGridSize(),pltwnd->z,100,NOFLOW);
		str_str=new strmr_strct(lapls);
		transfer_storage= new int[pltwnd->getGridSize()*pltwnd->getGridSize()];
		if (lapls->cudaStatus != cudaSuccess) {
			printf("cuda error!");
		}
		lapls->cpySlice(pltwnd->arr);
		str_str->cpySlice(transfer_storage);
		for(int i=0;i<pltwnd2->getGridSize();i++){
			for(int j=0;j<pltwnd2->getGridSize();j++){
				pltwnd2->arr[i+j*pltwnd->getGridSize()]=transfer_storage[i+j*pltwnd->getGridSize()];
			}
			//  std::cout<<'\n';
		}
		Fl::lock();
		fl_create_thread(fl_thread_calc, (Fl_Thread_Func*)thread_calc, lapls);
		pltwnd->show();
		pltwnd2->show();
		Fl::unlock();
		//lapls->cpySlice(pltwnd->arr);
		/// ((Fl_Window*)  pltwnd)->redraw();	  
	}
}
//Ёта функци€ срабатывает при нажатии кнопки старт, провер€ет состо€ние расчета, создает поток рисовани€ и поток расчета

void stop_cb(Fl_Widget *w,void* data){	
	Control_Window* tmp=(Control_Window*)data;
	if(tmp->calc_state & 1){
		((Fl_Window*)pltwnd)->hide();
		delete pltwnd;
		tmp->calc_state^=1;
		if(tmp->aftertasks & 1){
			tmp->aftertasks^=1;
			tmp->iin_cb( tmp->iin, data);
		}
	}
	if(tmp->aftertasks & 2){
		delete pltwnd;
		delete lapls;
		tmp->aftertasks ^=2;
	}
} 
//кнопка стоп и лучше eЄ не трогать

void thread_calc(void* p){
	cudaProfilerStart();
	double eps=.00001;
	double relEps=.00001;
	char* niv_check=new char[ctrlwnd->getGridSize()];
	char* relNiv_check=new char[ctrlwnd->getGridSize()];
	int nCube=ctrlwnd->getGridSize()*ctrlwnd->getGridSize()*ctrlwnd->getGridSize();
	for(int i=0;i<ctrlwnd->getGridSize();i++)niv_check[i]=0;
	for(int i=0;i<ctrlwnd->getGridSize();i++)relNiv_check[i]=0;
	int iterCounter=250000;
	output=new double[nCube];
	printf("start\n");
	//
	float intervalTime=0;

	//FILE* pF=fopen("field","w+");
	////for(int i=0;i<401*401*401;i++)output[i]=0;
	//for(;iterCounter<80000;iterCounter++)lapls->iteration(str_str,niv_check,eps,&intervalTime);
	//cudaMemcpy(output,lapls->dev_fi,sizeof(double)*nCube,cudaMemcpyDeviceToHost);

	//std::stringstream fileNameStream;
	//fileNameStream << "C:/Users/student1/Desktop/FOLDER_0/New folder/fiArr"<<ctrlwnd->getGridSize()<<"size" <<iterCounter<<'i'<<eps <<"eps.dat";
	//std::string fileName =fileNameStream.str();
	loadArray(output,nCube, "C:/Users/student1/Desktop/FOLDER_0/New folder/outEps8");

	//fileName.clear();
	//fileNameStream.clear();
	////double* outputCheck=new double[nCube];
	////loadArray(outputCheck,nCube, fileName);
	printf("conv =%d\n",niv_check[0]);
	printf("relConv =%d\n",relNiv_check[0]);
	printf("iter =%d\n",iterCounter);


	cudaMemcpy(lapls->dev_fi,output,sizeof(double)*ctrlwnd->getGridSize()*ctrlwnd->getGridSize()*ctrlwnd->getGridSize(),
		cudaMemcpyHostToDevice);
	cudaMemcpy(lapls->dev_fi_old,output,sizeof(double)*ctrlwnd->getGridSize()*ctrlwnd->getGridSize()*ctrlwnd->getGridSize(),
		cudaMemcpyHostToDevice);


	double k=0;


	k=( output[  
		ctrlwnd->getGridSize()/2
			+ctrlwnd->getGridSize()  *(ctrlwnd->getGridSize()  /2)
			+ctrlwnd->getGridSize()  *ctrlwnd->getGridSize()  *(ctrlwnd->getGridSize()  /2)]
	-output[
		ctrlwnd->getGridSize()/2
			+ctrlwnd->getGridSize()  *(ctrlwnd->getGridSize()/2)
			+ctrlwnd->getGridSize()  *ctrlwnd->getGridSize()  *(1  +ctrlwnd->getGridSize()  /2)]
	)/1 *450;

	str_str->count_report(/*fstream_h_N,*/k);
	lapls->cpySlice(pltwnd->arr);
	str_str->cpySlice(transfer_storage);

	delete[] output;
	float globalLapTime=0;
	int iterstep=100;
	relEps=eps;



	for(int i=0;i<ctrlwnd->getGridSize();i++)niv_check[0]*=niv_check[i];
	for(int i=0;i<ctrlwnd->getGridSize();i++)relNiv_check[0]*=relNiv_check[i];
	while(1){

		printf("conv =%d\n",niv_check[0]);
		printf("relConv =%d\n",relNiv_check[0]);
		printf("iter =%d\n",iterCounter);
		printf("%f\n",intervalTime);
		if(/*!niv_check[0]   &&*/ (ctrlwnd->calc_state==1)){


			for(int i=0;i<iterstep;i++,iterCounter++)lapls->iteration(str_str,niv_check,eps,&intervalTime);
			globalLapTime+=intervalTime;
			lapls->convergence(niv_check,eps);
			lapls->convergence(relNiv_check,relEps);

			for(int i=0;i<ctrlwnd->getGridSize();i++)niv_check[0]*=niv_check[i];
			for(int i=0;i<ctrlwnd->getGridSize();i++)relNiv_check[0]*=relNiv_check[i];

			if(niv_check[0]&&relNiv_check[0]){
				output=new double[nCube];

				//fileNameStream << "C:/Users/student1/Desktop/FOLDER_0/New folder/fiArr"<<ctrlwnd->getGridSize()<<"size" <<iterCounter/1000<<'i'<<eps <<"eps.dat";
				cudaMemcpy(output,lapls->dev_fi,sizeof(double)*nCube,cudaMemcpyDeviceToHost);
				//fileName =fileNameStream.str();

				//saveArray(output,nCube, "outEps3");

				//fileNameStream.clear();
				//fileName.clear();
				delete[] output;
				double* graphFi= new double[ctrlwnd->getGridSize()];
				lapls->GetFiCentral(graphFi);
				graphFi[0]=1;
				std::ofstream FiDat;
				FiDat.open("FiCentral.Dat",std::ofstream::out);
				FiDat.precision(16);
				for(int i=0;i<ctrlwnd->getGridSize();i++){
					FiDat<<graphFi[i]<<'\n';
				}
				FiDat.close();
				//eps *=0.1;
				//relEps=eps;
				/*lapls->MaxSearch(str_str->_states);
				str_str->cu_iterate(lapls->dev_fi);*/
				str_str->count_report(/*fstream_h_N,*/k);
			}
			lapls->cpySlice(pltwnd->arr);
			str_str->cpySlice(transfer_storage);


			//Fl::lock();
			//Fl::unlock();
			Fl::awake(pltwnd);
		}
		else {
			cudaProfilerStop();
			cudaDeviceReset();
			/*exit(0);*/
			return;
		}
	}
}
///это основна€ функци€ работ€юща€ в потоке расчета

#else
#  include <FL/fl_ask.H>

int main() {
	fl_alert("Sorry, threading not supported on this platform!");
}
#endif