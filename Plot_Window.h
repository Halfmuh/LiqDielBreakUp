#pragma once
#include <FL\Fl_Window.h>
#include "Control_Window.h"
#include <FL\fl_draw.H>
class Plot_Window :
	public Fl_Window
{
	int _grSz;
public:
	double* arr;
	int show_src;
	int z; 
	bool rotate;
	Control_Window* ctrl;
	Plot_Window(int W,int H,Control_Window *F,const char *L);
	static void plt_cb(Fl_Widget *w, void* data);
	///cpyfromdevice();
	int getGridSize();
	void draw() override;
};

