#pragma once
#include <string>
#include <sstream>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Int_Input.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Button.H>


class Control_Window :
	public Fl_Window
{
protected:
		int _grSz;
public:
	Fl_Int_Input *iin;
	Fl_Output *fout;
	Fl_Button *start;
	Fl_Button *stop;
	int calc_state;
	int aftertasks;
	inline void setGridSize();
	Control_Window(int W,int H,void(*sta)(Fl_Widget*,void*),void(*sto)(Fl_Widget*,void*),const char *L =0);
	static void ctrl_cb(Fl_Widget* w,void* data);
	static void iin_cb(Fl_Widget* w,void* data);
	int getGridSize();
	~Control_Window(void);
};

