#include <iostream>
#include "config_loader.h"
#include "sample.h"
#include "mass.h"

using namespace std;
using namespace Eigen;

void onlineSim()
{
    cout << "onlineSim" << endl;
    VectorXd initPose, initVel;
    Utility::setStateAt(0, initPose, initVel);
    cout << "initial pose:" << endl;
    cout << initPose << endl;
    cout << "initial vel:" << endl;
    cout << initVel << endl;
    vector<shared_ptr<Sample>> samples{ make_shared<Sample>(initPose, initVel) };
    for (int i = 0; i < walk.size(); ++i)
    {
	cout << "i = " << i << endl;
	cout << "frag " << walk[i] << endl;
	ControlFragment &frag = frags[walk[i]];
	shared_ptr<Sample> prev = samples.back();
	VectorXd action = frag.m * makeState(prev->resultPose, prev->resultVel) + frag.a;
	VectorXd offset = action2offset(action);
	bool useID = true;
	shared_ptr<Sample> sample = make_shared<Sample>(prev, frag, offset, action, simulators[omp_get_thread_num()], useID);
	samples.push_back(sample);
	cout << sample->cost << endl;
	if (sample->cost > Config::failThreshold)
	    break;
    }
    ofstream output(Config::outputFileName + "_online.txt");
    for (const Eigen::VectorXd &v: Utility::bvhs[omp_get_thread_num()].frameToEulerAngle(samples.back()->getTrajectory())) // used bvh4window.frame instead?
	output << v.transpose() << std::endl;
    output.close();
}
