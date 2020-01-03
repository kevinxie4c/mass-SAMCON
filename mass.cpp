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
vector<VectorXd> forces;
Timer timer;


int main(int argc, char* argv[])
{
#ifdef NDEBUG
    cout << "mass (release)" << endl;
#else
    cout << "mass (debug)" << endl;
#endif
    cout << "use inverse dynamics" << endl;
    cout << "notImprove: use accumulative error to decide" << endl;
    if (argc > 1)
    {
	if (argc > 4)
	{
	    cout << argv[0] << " [task_file] [use_mass] [rounds]" << endl;
	    exit(0);
	}
	if (argc >= 3)
	    useMass = stoi(argv[2]);
	if (argc >= 4)
	    rounds = stoi(argv[3]);
	taskFileName = argv[1];
    }
    cout << "useMass = " << useMass << endl;
    Config::load(taskFileName);
    Utility::init();
    if (Config::endFrame > Utility::mbvh.frame.size())
	Config::endFrame = Utility::mbvh.frame.size();
    vector<VectorXd> frame(Utility::mbvh.frame.cbegin() + Config::startFrame, Utility::mbvh.frame.cbegin()+ Config::endFrame);
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

    if (Config::useInverseForce)
    {
	forces = Utility::readVectorXdListFrom(Config::forceFileName);
	forces = std::vector<Eigen::VectorXd>(forces.cbegin() + Config::startFrame, forces.cbegin() + Config::endFrame);
    }
    else
	forces = std::vector<Eigen::VectorXd>(Config::endFrame - Config::startFrame, Eigen::VectorXd::Zero(Utility::ndof));
    setUpFrags(useMass);

    for (size_t i = 0; i < Config::loopNum; ++i)
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

    cout << "duration: " << timer.durationToString() << endl;;
    if (Config::generateSamplesFile)
    {
	if (!Utility::dirExists("samples"))
	    Utility::createDir("samples");
    }
    while (rounds--)
    {
	Timer t;
	refine(useMass);
	cout << "round duration: " << t.durationToString() << endl;
    }

    return 0;
}
