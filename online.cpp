#include <iostream>
#include "config_loader.h"
#include "utility.h"
#include "sample.h"
#include "mass.h"

using namespace std;
using namespace Eigen;

void onlineSim()
{
    cout << "onlineSim" << endl;
    VectorXd initPose, initVel;
    if (Config::initStateFileName == "")
	Utility::setStateAt(0, initPose, initVel);
    else
    {
	if (Utility::fileGood(Config::initStateFileName))
	{
	    vector<VectorXd> list = Utility::readVectorXdListFrom(Config::initStateFileName);
	    initPose = list[0];
	    initVel = list[1];
	}
	else
	{
	    cout << "cannot open file " << Config::initStateFileName << endl;
	    exit(0);
	}
    }
    cout << "initial pose:" << endl;
    cout << initPose << endl;
    cout << "initial vel:" << endl;
    cout << initVel << endl;
    vector<shared_ptr<Sample>> samples{ make_shared<Sample>(initPose, initVel) };
    Simulator &sim = simulators[omp_get_thread_num()];
    sim.setPose(initPose, initVel);
    SkeletonPtr skeleton = sim.skeleton;
    for (int i = 0; i < walk.size(); ++i)
    {
	cout << "i = " << i << endl;
	cout << "frag " << walk[i] << endl;
	ControlFragment &frag = frags[walk[i]];
	/*
	shared_ptr<Sample> prev = samples.back();
	VectorXd action = frag.m * makeState(prev->resultPose, prev->resultVel) + frag.a;
	VectorXd offset = action2offset(action);
	bool useID = true;
	shared_ptr<Sample> sample = make_shared<Sample>(prev, frag, offset, action, simulators[omp_get_thread_num()], useID);
	samples.push_back(sample);
	cout << sample->cost << endl;
	if (sample->cost > Config::failThreshold)
	    break;
	*/

	// TODO: rewrite according to sample.cpp and simulator.cpp (e.g., to deal with option such as stablePD).
	VectorXd action = frag.m * makeState(sim.skeleton->getPositions(), sim.skeleton->getVelocities()) + frag.a;
	VectorXd offset = action2offset(action);
	for (size_t j = 0; j < 6; ++j)
	    offset[j] = 0;
	VectorXd ref = frag.tracked + offset;
	for (size_t j = 0; j < Config::groupNum; ++j)
	{
	    // Stable PD
	    Eigen::VectorXd q = skeleton->getPositions();
	    Eigen::VectorXd dq = skeleton->getVelocities();
	    Eigen::MatrixXd invM = (skeleton->getMassMatrix() + Utility::mKd * skeleton->getTimeStep()).inverse();
	    Eigen::VectorXd p = -Utility::mKp * skeleton->getPositionDifferences(q + dq * skeleton->getTimeStep(), ref);
	    Eigen::VectorXd d = -Utility::mKd * dq;
	    Eigen::VectorXd qddot = invM * (-skeleton->getCoriolisAndGravityForces() + p + d + skeleton->getConstraintForces());
	    Eigen::VectorXd force = p + d -Utility::mKd * qddot * skeleton->getTimeStep();
	    force += frag.iforces[j];
	    //Eigen::VectorXd force = iforce[i];
	    forces.push_back(force);
	    for (size_t k = 0; k < Config::stepPerFrame; ++k)
	    {
		skeleton->setForces(force);
		sim.world->step();
	    }
	}
	Vector3d zmp;
	Utility::ErrorTerms et;
	double cost = Utility::costFunc(sim.skeleton, frag, zmp, et);
	cout << cost << endl;
	if (cost > Config::failThreshold)
	    break;

    }
    ofstream output(Config::outputFileName + "_online.txt");
    for (const Eigen::VectorXd &v: Utility::bvhs[omp_get_thread_num()].frameToEulerAngle(samples.back()->getTrajectory())) // used bvh4window.frame instead?
	output << v.transpose() << std::endl;
    output.close();
}
