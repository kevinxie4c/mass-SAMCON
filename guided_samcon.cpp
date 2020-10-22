#include <iostream>
#include <queue>
#include <omp.h>
#include "config_loader.h"
#include "eigenmvn.h"
#include "sample.h"
#include "MyWindow.h"
#include "mass.h"

using namespace std;
using namespace Eigen;


VectorXd action2offset(VectorXd &action)
{
    assert(action.rows() == 11);
    Vector3d waist = action.segment(0, 3);
    Vector3d leftHip = action.segment(3, 3);
    Vector3d rightHip = action.segment(6, 3);
    VectorXd leftKnee = action.row(9);
    VectorXd rightKnee = action.row(10);
    SkeletonPtr skeleton = Utility::bvhs[omp_get_thread_num()].skeleton;
    skeleton->setPositions(Eigen::VectorXd::Zero(Utility::ndof));
    Joint* joint = skeleton->getJoint(Utility::waistIndex);
    joint->setPositions(waist);
    joint = skeleton->getJoint(Utility::leftHipIndex);
    joint->setPositions(leftHip);
    joint = skeleton->getJoint(Utility::rightHipIndex);
    joint->setPositions(rightHip);
    joint = skeleton->getJoint(Utility::leftKneeIndex);
    joint->setPositions(leftKnee);
    joint = skeleton->getJoint(Utility::rightKneeIndex);
    joint->setPositions(rightKnee);
    return skeleton->getPositions();
}

static const double lambda = 1e-6;

// re-calculate frags
void recalcFrag(const vector<vector<VectorXd>> &statesOf, const vector<vector<VectorXd>> &actionsOf)
{
    for (size_t i = 0; i < frags.size(); ++i)
    {
	const vector<VectorXd> &sv = statesOf[i];
	MatrixXd sm(18, sv.size());
	for (size_t j = 0; j < sv.size(); ++j)
	    sm.col(j) = sv[j];
	VectorXd smean = sm.rowwise().mean();
	sm.colwise() -= smean;
	sm.transposeInPlace();
	
	const vector<VectorXd> &av = actionsOf[i];
	MatrixXd am(11, av.size());
	for (size_t j = 0; j < av.size(); ++j)
	    am.col(j) = av[j];
	VectorXd amean = am.rowwise().mean();
	am.colwise() -= amean;
	am.transposeInPlace();

	assert(sv.size() == av.size());

	frags[i].m = ((sm.transpose() * sm + lambda * MatrixXd::Identity(18, 18)).inverse() * (sm.transpose() * am)).transpose();
	frags[i].a = amean - frags[i].m * smean;

	MatrixXd tm = am - sm * frags[i].m.transpose();
	// 1.0 / sv.size() or 1.0 / walk.size() ?
	frags[i].sigma = ((1.0 / sv.size()) * (tm.transpose() * tm).diagonal()).asDiagonal();
    }
}

/*
 *        f[i]
 * s[i] -------> s[i+1]
 *
 */

void guidedSAMCON()
{
    static size_t counter = 0;
    const static size_t initRetreat = 10;
    const static size_t maxRetreat = 50;
    ++counter;
    cout << "guidedSAMCON " << counter << endl;
    VectorXd initPose, initVel;
    Utility::setStateAt(0, initPose, initVel);
    cout << "initial pose:" << endl;
    cout << initPose << endl;
    cout << "initial vel:" << endl;
    cout << initVel << endl;
    vector<vector<shared_ptr<Sample>>> savedSamples;
    vector<shared_ptr<Sample>> tmpList;
    for (size_t i = 0; i < Config::saveNum; ++i)
	tmpList.push_back(make_shared<Sample>(initPose, initVel));
    savedSamples.push_back(tmpList);
    vector<size_t> retreats(walk.size(), initRetreat);
    vector<EigenMultivariateNormal<double>> normalDists;
    for (ControlFragment &cf: frags)
    {
	//cout << cf.sigma << endl;
	normalDists.push_back(EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(11), cf.sigma));
    }
    for (int i = 0; i < walk.size(); ++i)
    {
	cout << "i = " << i << endl;
	cout << "frag " << walk[i] << endl;
	ControlFragment &frag = frags[walk[i]];
	EigenMultivariateNormal<double> &dist = normalDists[walk[i]];
	// we use savedSamples[i - 1 + 1]
	// save sample at savedSamples[i + 1]
	vector<shared_ptr<Sample>> &samples = savedSamples[i];
	priority_queue<shared_ptr<Sample>, vector<shared_ptr<Sample>>, CostCmp> queue;
#pragma omp parallel for
	for (size_t j = 0; j < samples.size(); ++j)
	{
	    shared_ptr<Sample> &sample = samples[j];
	    VectorXd mean = frag.m * makeState(sample->resultPose, sample->resultVel) + frag.a;
	    //cout << "mean:" << endl;
	    //cout << mean << endl;
	    for (size_t k = 0; k < Config::sampleNum / samples.size(); ++k)
	    {
		VectorXd action = mean + dist.samples(1);
		//cout << "action:" << endl;
		//cout << action << endl;
		VectorXd offset = action2offset(action);
		//cout << "offset:" << endl;
		//cout << offset << endl;
		bool useID = true;
		shared_ptr<Sample> t = make_shared<Sample>(sample, frag, offset, action, simulators[omp_get_thread_num()], useID);
#pragma omp critical (queue_section)
		{
		    queue.push(t);
		    if (queue.size() > Config::saveNum)
			queue.pop();
		}
	    }
	}
	vector<shared_ptr<Sample>> tmp;
	while (!queue.empty())
	{
	    const shared_ptr<Sample> &t = queue.top();
	    if (!isnan(t->cost))
		tmp.push_back(queue.top());
	    queue.pop();
	}
	cout << tmp.back()->cost << endl;
	if (tmp.back()->cost > Config::failThreshold)
	{
	    cout << "fail" << endl;
	    std::shared_ptr<Sample> minSample = nullptr;
	    double min = DBL_MAX;
	    // last one is savedSamples[i_end - 1 + 1]
	    for (std::shared_ptr<Sample> &sample: savedSamples[i])
	    {
		if (sample->totalCost < min)
		{
		    min = sample->totalCost;
		    minSample = sample;
		}
	    }
	    vector<shared_ptr<const Sample>> minSamplesList = minSample->getMinSamplesList();
	    if (Config::showWindow)
	    {
		MyWindow::bvh4window.frame.clear();
		for (std::shared_ptr<const Sample> s: minSamplesList)
		{
		    for (auto it: s->trajectory)
			MyWindow::bvh4window.frame.push_back(it);
		}
	    }

	    size_t t = i;
	    i -= retreats[i];
	    retreats[t] *= 2;
	    if (retreats[t] > maxRetreat)
		retreats[t] = maxRetreat;
	    if (i < 0) i = 0;
	    i -= 1;
	    continue;
	}
	if (retreats[i] > initRetreat)
	    --retreats[i];
	assert(i + 1 >= 0);
	if ((size_t)(i + 1) < savedSamples.size()) // save at i + 1 instead of i
	    savedSamples[i + 1] = tmp;
	else
	    savedSamples.push_back(tmp);
    }
    std::shared_ptr<Sample> minSample = nullptr;
    double min = DBL_MAX;
    // last one is savedSamples[i_end - 1 + 1]
    for (std::shared_ptr<Sample> &sample: savedSamples.back())
    {
	if (sample->totalCost < min)
	{
	    min = sample->totalCost;
	    minSample = sample;
	}
    }
    std::cout << "trial: " << counter << std::endl;
    vector<shared_ptr<const Sample>> minSamplesList = minSample->getMinSamplesList();
    vector<vector<VectorXd>> statesOf(frags.size()), actionsOf(frags.size());
    for (size_t i = 0; i < walk.size(); ++i)
    {
	size_t index = walk[i];
	statesOf[index].push_back(makeState(minSamplesList[i]->resultPose, minSamplesList[i]->resultVel));
	actionsOf[index].push_back(minSamplesList[i + 1]->kernel);
    }
    recalcFrag(statesOf, actionsOf);
    std::ofstream output;
    for (size_t i = 0; i < frags.size(); ++i)
    {
	ControlFragment &frag = frags[i];

	output.open("out_frags_" + std::to_string(i) + "_m.txt");
	output << frag.m << endl;
	output.close();

	output.open("out_frags_" + std::to_string(i) + "_a.txt");
	output << frag.a << endl;
	output.close();

	output.open("out_frags_" + std::to_string(i) + "_sigma.txt");
	output << frag.sigma << endl;
	output.close();
    }

    // compute the average position and velocity (for setting up initial position and velocity)
    vector<vector<VectorXd>> positionOf(frags.size()), velocityOf(frags.size());
    for (size_t i = 0; i < walk.size(); ++i)
    {
	size_t index = walk[i];
	positionOf[index].push_back(minSamplesList[i]->resultPose);
	velocityOf[index].push_back(minSamplesList[i]->resultVel);
    }
    for (size_t i = 0; i < positionOf.size(); ++i)
    {
	// output initPose, initVel for frags[i]
	VectorXd p = VectorXd::Zero(initPose.rows());
	VectorXd v = VectorXd::Zero(initVel.rows());
	for (const VectorXd &u: positionOf[i])
	    p += u;
	p /= positionOf[i].size();
	for (const VectorXd &u: velocityOf[i])
	    v += u;
	v /= velocityOf[i].size();
	output.open("init_state_" + to_string(i) + ".txt");
	output << p.transpose() << endl;
	output << v.transpose() << endl;
	output.close();
    }

    if (Config::showWindow)
    {
	MyWindow::bvh4window.frame.clear();
	for (std::shared_ptr<const Sample> s: minSamplesList)
	{
	    for (auto it: s->trajectory)
		MyWindow::bvh4window.frame.push_back(it);
	}
    }
    output.open(Config::outputFileName + "_guided_"  + std::to_string(counter) + ".txt");
    for (const Eigen::VectorXd &v: Utility::bvhs[omp_get_thread_num()].frameToEulerAngle(minSample->getTrajectory())) // used bvh4window.frame instead?
	output << v.transpose() << std::endl;
    output.close();
#ifndef NDEBUG
    vector<VectorXd> com, mmt;
    for (std::shared_ptr<const Sample> s: minSamplesList)
    {
	if (!Config::onlyLogAndFinal)
	{
	    for (auto it: s->com)
		com.push_back(it);
	    for (auto it: s->mmt)
		mmt.push_back(it);
	}
    }
    if (!Config::onlyLogAndFinal)
    {
	output.open("com" + std::to_string(counter) + ".txt");
	for (const Eigen::VectorXd &v: com)
	    output << v.transpose() << std::endl;
	output.close();
	output.open("mmt" + std::to_string(counter) + ".txt");
	for (const Eigen::VectorXd &v: mmt)
	    output << v.transpose() << std::endl;
	output.close();
    }
#endif

}
