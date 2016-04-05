#include "Control_window.h"


Control_Window::Control_Window(int W,int H,void(*sta)(Fl_Widget*,void*), void(*sto)(Fl_Widget*,void*),const char *L):Fl_Window(W,H,L){
	  calc_state=0;
	  aftertasks=0;
	  _grSz=301;
	  begin();
	  callback(ctrl_cb,NULL);
	  //iin = new Fl_Int_Input(100,10,100,30,"Grid side");
	  //iin->value("20");
	  //fout= new Fl_Output(100,50,100,30,"Grid side");
	  //fout->value("20");
	  start = new Fl_Button(100,90,50,20,"start");
	  start->callback(sta,this);
	  stop = new Fl_Button(100,120,50,20,"stop");
	  stop->callback(sto,this);
	  //iin->callback(iin_cb, this);
	  //iin->when(FL_WHEN_RELEASE);
	  //fout->when(FL_WHEN_CHANGED);
	  end();
}

void Control_Window::iin_cb(Fl_Widget *w, void* data){
	Control_Window* tmp=(Control_Window*)data;
	tmp->setGridSize();
	
}

inline void Control_Window::setGridSize(){
	if (calc_state & 1) aftertasks |= 1;
	//else {
	//	char a[10];
	//	for(int i=0;i<10;i++)a[i]=0;
	//	std::stringstream stream;
	//	stream.str(iin->value());
	//	stream>>_grSz;
	//	stream<<_grSz;
	//	stream>>a;
	//	fout->value(a);                
	// }
}

void Control_Window::ctrl_cb(Fl_Widget* w,void* data){
	exit(0);
}

int Control_Window::getGridSize(){
return _grSz;
}
Control_Window::~Control_Window(void)
{
}
