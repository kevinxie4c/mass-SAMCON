#include <iostream>
#include <pthread.h>
#include "config_loader.h"
#include "utility.h"
#include "simulator.h"
#include "control_fragment.h"
#include "BVHData.h"
#include "cmaes.h"
#include "MyWindow.h"
#include "mass.h"

using namespace std;
using namespace Eigen;

size_t numFrag;
string taskFileName("task_walk.txt");
size_t rounds = 5;
bool useMass = true;
vector<size_t> walk;
vector<ControlFragment> frags;
vector<VectorXd> initMean;
vector<Simulator> simulators;


int main(int argc, char* argv[])
{
#ifdef NDEBUG
    cout << "mass (release)" << endl;
#else
    cout << "mass (debug)" << endl;
#endif
    if (argc > 1)
    {
	if (argc > 4)
	{
	    cout << argv[0] << " [bvh_file] [use_mass] [rounds]" << endl;
	    exit(0);
	}
	if (argc >= 3)
	    useMass = stoi(argv[2]);
	if (argc >= 4)
	    rounds = stoi(argv[3]);
	taskFileName = argv[1];
    }
    Config::load(taskFileName);
    Utility::init();
    vector<VectorXd> frame(Utility::mbvh.frame.cbegin() + Config::startFrame, Utility::mbvh.frame.cend());
    for (BVHData &bvh: Utility::bvhs)
	bvh.frame = frame;
    Utility::mbvh.frame = frame;
    cout << "frag frame size: " << frame.size() << endl;
    numFrag = frame.size() / Config::groupNum;
    cout << "numFrag: " << numFrag << endl;

    pthread_t thread;
    MyWindow window;
    WorldPtr world(new World);
    if (Config::showWindow)
    {
	world->addSkeleton(MyWindow::bvh4window.skeleton);
	world->setTimeStep(0.1);
	window.setWorld(world);
	glutInit(&argc, argv);
	window.initWindow(640, 480, "test");
	std::cout << "q: frame - 1" << std::endl;
	std::cout << "w: frame + 1" << std::endl;
	std::cout << "a: frame - 10" << std::endl;
	std::cout << "s: frame - 10" << std::endl;
	std::cout << "p: play" << std::endl;
	pthread_create(&thread, NULL, (void* (*)(void*))glutMainLoop, NULL);
    }

    setUpFrags(useMass);

    for (size_t j = 0; j < numFrag; ++j)
	walk.push_back(j);

#ifdef _OPENMP
    cout << "_OPENMP " << _OPENMP << endl;
#else
    cout << "NO OPENMP" << endl;
#endif
    cout << "omp_get_max_threads() = " << omp_get_max_threads() << endl;
    for (size_t i = 0; i < omp_get_max_threads(); ++i)
	simulators.push_back(Simulator(Utility::bvhs[omp_get_thread_num()]));

    // full rank if we are using original algorithm
    if (!useMass)
	Config::rank = Utility::ndof;
    if (initMean.empty())
	for (size_t i = 0; i < walk.size(); ++i)
	    initMean.push_back(Eigen::VectorXd::Zero(Config::rank));

    while (rounds--)
	refine(useMass);

    return 0;
}
