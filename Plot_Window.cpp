#include "Plot_Window.h"
#include <stdio.h>
#include <iostream>
///constructor
Plot_Window::Plot_Window(int W,int H,Control_Window *F,const char *L=0):Fl_Window(W,H,L){
	ctrl=F;
	//
	_grSz=F->getGridSize();
	///
	z= (F->getGridSize())/2;
	rotate=0;
	show_src=0;
	arr =new double[ _grSz * _grSz ];
///	fl_create_thread(calc_thread, (Fl_Thread_Func *)calculation, pltwnd);
}


////// work to do:   make this CallBack beautiful!
void Plot_Window::plt_cb(Fl_Widget *w, void* data){
	Control_Window* tmp=(((Plot_Window*)data)->ctrl);
	((Fl_Window*)data)->hide();
	tmp->calc_state^=1;
	if(tmp->aftertasks & 1){
		tmp->aftertasks^=1;
		tmp->iin_cb( tmp->iin, data);
	}
	tmp->aftertasks ^=2;
  }
/// razuznat' razmer setki
int Plot_Window::getGridSize(){
	return _grSz;
}
///risowashki
void Plot_Window::draw(){
	Fl_Color col;
	int sz=getGridSize();
	static int n=0;
	int x_max=w();int y_max=h();
	int dy=y_max/sz ;int dx=x_max/sz ;
	  int dx2=(dx>>1);int dy2=(dy>>1);
	  int y=0;int j;
	  int x=-dx2;int i;

	  for(i=0;i<sz;x+=dx,i++){
		  for(j=0,y=-dy2;j<sz;y+=dy,j++){
			/*	std::cout<<arr[i+j*sz]<<" ";*/
			if (arr[i+j*sz]<=1){
				if(arr[i+j*sz]<=0.25)col=fl_color_average(FL_BLUE,FL_BLACK,(float)arr[i+j*sz]*4);
				else if(arr[i+j*sz]<=0.5)col=fl_color_average(FL_GREEN,FL_BLUE,(float)arr[i+j*sz]*4 -1);
				else if(arr[i+j*sz]<=0.75)col=fl_color_average(FL_YELLOW,FL_GREEN,(float)arr[i+j*sz]*4-2);
				else col=fl_color_average(FL_RED,FL_YELLOW,(float)arr[i+j*sz]*4-3);
			}
			else if(arr[i+j*sz]<30 && arr[i+j*sz]>1)col=FL_RED;
			else if(arr[i+j*sz]<60 && arr[i+j*sz]>31)col=fl_color_average(FL_BLACK,FL_RED,(float)0.5);
			else if(arr[i+j*sz]==101) col =FL_BLACK;
			else col=FL_GREEN;
			fl_rectf(x,y,dx,dy,col);
		  }
		  /*std::cout<<'\n';*/
	  }
	  ///printf("draw");
}

