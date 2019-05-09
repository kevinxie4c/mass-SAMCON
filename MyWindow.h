#ifndef MYWINDOW_H
#define MYWINDOW_H

#include <dart/gui/gui.hpp>
#include "BVHData.h"

class MyWindow: public dart::gui::SimWindow
{
    public:
	bool toggle = false;
	int index = 0;

	static BVHData bvh4window;

	MyWindow();

	void keyboard(unsigned char key, int x, int y) override;
	void timeStepping() override;
	
};

#endif
