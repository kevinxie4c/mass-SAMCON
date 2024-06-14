#include "config_loader.h"
#include "MyWindow.h"

BVHData MyWindow::bvh4window;

MyWindow::MyWindow()
{
    bvh4window.loadBVH(Config::bvhFileName, Config::geometryConfigFileName, Config::hingeJointListFileName, Config::scale);
}

void MyWindow::keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
	case 'q':
	    --index;
	    break;
	case 'w':
	    ++index;
	    break;
	case 'a':
	    index -= 10;
	    break;
	case 's':
	    index += 10;
	    break;
	case 'p':
	    toggle = !toggle;
	    break;
	default:
	    SimWindow::keyboard(key, x, y);
    }
    while (index < 0)
	index += bvh4window.frame.size();
    index = index % bvh4window.frame.size();
    bvh4window.setPositionAt(index);
    //bvhRef.setPositionAt(index);
}

void MyWindow::timeStepping()
{
    if (toggle)
    {
	++index;
	index = index % bvh4window.frame.size();
	bvh4window.setPositionAt(index);
	//bvhRef.setPositionAt(index);
    }
}
