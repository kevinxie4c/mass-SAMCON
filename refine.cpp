#include <iostream>
#include <queue>
#include <omp.h>
#include "config_loader.h"
#include "cmaes.h"
#include "sample.h"
#include "MyWindow.h"
#include "mass.h"

using namespace std;
using namespace Eigen;

/*
 *        f[i]
 * s[i] -------> s[i+1]
 *
 */

void refine(bool useMass)
{
    static size_t counter = 0;
    ++counter;
    cout << "refine " << counter << endl;
    static double sigmaMax = Config::initSigma;
    size_t trial = 0;
    vector<WeirdCMAES> cmaes;
    cout << initMean[0] << endl;
    for (size_t i = 0; i < walk.size(); ++i)
	cmaes.push_back(WeirdCMAES(Config::rank, Config::sampleNum, Config::saveNum, Config::initSigma, initMean[i]));
    vector<size_t> generation(walk.size(), 0);
    vector<size_t> notImprove(walk.size(), 0);
    vector<double> minCost(walk.size(), DBL_MAX);
    vector<size_t> maxHeight(walk.size(), 0);
    vector<double> minAccCost(walk.size(), DBL_MAX);
    vector<vector<shared_ptr<Sample>>> savedSamples, backupSamples;
    double backupMin = DBL_MAX;
    VectorXd initPose, initVel;
    // assume that we start from first frame
    Utility::setStateAt(0, initPose, initVel);
    cout << "initial pose:" << endl;
    cout << initPose << endl;
    cout << "initial vel:" << endl;
    cout << initVel << endl;
    vector<shared_ptr<Sample>> tmpList;
    for (size_t i = 0; i < Config::saveNum; ++i)
	tmpList.push_back(make_shared<Sample>(initPose, initVel));
    savedSamples.push_back(tmpList);
    size_t i_begin = 0, i_end, i_begin_backup = 0;
    size_t trialTimes = 0;
    while (i_begin < walk.size())
    {
	i_end = i_begin + Config::slidingWindow;
	if (i_end > walk.size())
	    i_end = walk.size();
	bool failed = false;
	bool failedAtBegin = false;
	for (size_t i = i_begin; i < i_end; ++i)
	{
	    cout << "i = " << i << endl;
	    cout << "frag " << walk[i] << endl;
	    ControlFragment &frag = frags[walk[i]];
	    // we use savedSamples[i - 1 + 1]
	    // save sample at savedSamples[i + 1]
	    vector<shared_ptr<Sample>> &samples = savedSamples[i];
	    priority_queue<shared_ptr<Sample>, vector<shared_ptr<Sample>>, CostCmp> queue;
#pragma omp parallel for
	    for (size_t j = 0; j < samples.size(); ++j)
	    {
		shared_ptr<Sample> &sample = samples[j];
		for (size_t k = 0; k < Config::sampleNum / samples.size(); ++k)
		{
		    VectorXd kernel = cmaes[i].getSample();
		    VectorXd delta;
		    if (useMass)
			delta = frag.transformation.leftCols(Config::rank + 6).rightCols(Config::rank) * kernel;
		    else
			delta = kernel;
		    shared_ptr<Sample> t = make_shared<Sample>(sample, frag, delta, kernel, simulators[omp_get_thread_num()]);
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
		i_end = i;
		cout << "fail" << endl;
		failed = true;
		if (i == i_begin)
		    failedAtBegin = true;
		break;
	    }
	    if (tmp.back()->cost < minCost[i])
	    {
		minCost[i] = tmp.back()->cost;
		notImprove[i] = 0;
	    }
	    else
		++notImprove[i];
	    if (i + 1 < savedSamples.size()) // save at i + 1 instead of i
		savedSamples[i + 1] = tmp;
	    else
		savedSamples.push_back(tmp);
	}
	std::shared_ptr<Sample> minSample = nullptr;
	double min = DBL_MAX;
	// last one is savedSamples[i_end - 1 + 1]
	for (std::shared_ptr<Sample> &sample: savedSamples[i_end])
	{
	    if (sample->totalCost < min)
	    {
		min = sample->totalCost;
		minSample = sample;
	    }
	}
	++trial;
	cout << "duration: " << timer.durationToString() << endl;;
	std::cout << "trial: " << counter << " - " << trial << std::endl;
	std::vector<std::shared_ptr<const Sample>> minSamplesList = minSample->getMinSamplesList();
	std::vector<Eigen::VectorXd> minTrajectory;
	for (std::shared_ptr<const Sample> s: minSamplesList)
	{
	    std::cout << s->cost << " ";
	    minTrajectory.insert(minTrajectory.end(), s->trajectory.begin(), s->trajectory.end());
	}
	std::cout << endl;
	if (Config::showWindow)
	{
	    MyWindow::bvh4window.frame = minTrajectory;
	}
	if (!Config::onlyLogAndFinal)
	{
	    std::ofstream output;
	    if (!Config::onlyLogAndFinal)
		output.open(Config::outputFileName + std::to_string(counter) + "_" + std::to_string(trial) + ".txt");
	    for (const Eigen::VectorXd &v: Utility::bvhs[omp_get_thread_num()].frameToEulerAngle(minTrajectory)) // used bvh4window.frame instead?
		output << v.transpose() << std::endl;
	    output.close();
#ifndef NDEBUG
	    vector<VectorXd> com, mmt, forces;
	    for (std::shared_ptr<const Sample> s: minSamplesList)
	    {
		if (!Config::onlyLogAndFinal)
		{
		    for (auto it: s->com)
			com.push_back(it);
		    for (auto it: s->mmt)
			mmt.push_back(it);
		    for (auto it: s->forces)
			forces.push_back(it);
		}
	    }
	    output.open("com" + std::to_string(counter) + "_" + std::to_string(trial) + ".txt");
	    for (const Eigen::VectorXd &v: com)
		output << v.transpose() << std::endl;
	    output.close();
	    output.open("mmt" + std::to_string(counter) + "_" + std::to_string(trial) + ".txt");
	    for (const Eigen::VectorXd &v: mmt)
		output << v.transpose() << std::endl;
	    output.close();
	    output.open("forces" + std::to_string(counter) + "_" + std::to_string(trial) + ".txt");
	    for (const Eigen::VectorXd &v: forces)
		output << v.transpose() << std::endl;
	    output.close();
#endif
	}
	// current max step: i_end - 1
	// previous max step: backupSamples.size() - 1 - 1 (sample of step i is saved at samples[i + 1])
	// i_end - 1 > backupSamples.size() - 1 -1  ==>  i_end + 1 > backupSamples.size()
	if (i_end + 1 > backupSamples.size() || (i_end + 1 == backupSamples.size() && min < backupMin))
	{
	    for (size_t i = i_begin; i < i_begin + Config::updateWindow && i < i_end; ++i) // CMA-ES
	    {
		std::vector<std::shared_ptr<Sample>> &tmp = savedSamples[i + 1];
		std::priority_queue<std::shared_ptr<Sample>, std::vector<std::shared_ptr<Sample>>, SampleComparison> tQueue(tmp.cbegin(), tmp.cend());
		//Eigen::MatrixXd X(Config::rank + Config::dRank * (counter - 1), Config::saveNum);
		Eigen::MatrixXd X(Config::rank, Config::saveNum);
		size_t h = 0;
		double cost = 0;
		for (size_t j = 0; j < Config::saveNum; ++j)
		{
		    std::shared_ptr<Sample> s = tQueue.top();
		    h += s->height;
		    cost += s->accumCost;
		    X.col(j) = s->kernel;
		    tQueue.pop();
		}
		std::cout << i << " " << h << " " << cost << std::endl;
		//if (h > maxHeight[i] || (h == maxHeight[i] && cost < minAccCost[i])) // do CMA-ES only if the samples are better
		{
		    maxHeight[i] = h;
		    minAccCost[i] = cost;
		    cmaes[i].update(X);
		    if (cmaes[i].sigma > sigmaMax)
			cmaes[i].sigma = sigmaMax;
		}
	    }
	    backupSamples.clear();
	    for (size_t i = 0; i < i_end + 1; ++i)
		backupSamples.push_back(savedSamples[i]);
	    backupMin = min;
	}
	else
	{
	    /*
	    savedSamples.clear();
	    for (auto i: backupSamples)
		savedSamples.push_back(i);
	    */
	    savedSamples = backupSamples;
	}
	std::cout << std::endl;
	for (size_t i = 0; i < walk.size(); ++i)
	    std::cout << maxHeight[i] << " ";
	std::cout << std::endl;
	for (size_t i = 0; i < walk.size(); ++i)
	    std::cout << minAccCost[i] << " ";
	std::cout << std::endl;
	for (size_t i = 0; i < walk.size(); ++i)
	{
	    generation[i] = cmaes[i].counteval / cmaes[i].lambda;
	    std::cout << generation[i] << " ";
	}
	std::cout << std::endl;
	for (size_t i = 0; i < walk.size(); ++i)
	    std::cout << cmaes[i].sigma << " ";
	std::cout << std::endl;
	/*
	if (generation[i_begin] < Config::trialMin)
	{
	    if (++trialTimes > Config::trialMin)
	    {
		++i_begin;
		trialTimes = 0;
	    }
	    continue;
	}
	*/
	
	if (failedAtBegin)
	{
	    // rollback to the last successfull trial if failed at the beginning
	    i_begin = i_begin_backup;
	}
	else
	{
	    if (!failed)
		i_begin_backup = i_begin;
	    while (i_begin < walk.size() && (generation[i_begin] > Config::trialMax || notImprove[i_begin] > Config::notImproveMax || minCost[i_begin] < Config::goodEnough))
	    {
		++i_begin;
	    }
	}
	trialTimes = 0;
    }
    std::shared_ptr<Sample> minSample = nullptr;
    double min = DBL_MAX;
    for (std::shared_ptr<Sample> sample: savedSamples.back())
    {
	if (sample->totalCost < min)
	{
	    min = sample->totalCost;
	    minSample = sample;
	}
    }
    
    // noise reduction
    for (size_t i = 1; i < savedSamples.size(); ++i)
    {
	size_t h = 0;
	initMean[i - 1] = Eigen::VectorXd::Zero(Config::rank);
	Eigen::VectorXd mean = Eigen::VectorXd::Zero(Utility::ndof);
	for (const std::shared_ptr<Sample> &s: savedSamples[i])
	{
	    if (s->height >= 4)
	    {
		initMean[i - 1] += s->kernel * s->height;
		h += s->height;
		mean += s->delta * s->height;
	    }
	}
	if (h > 0)
	{
	    initMean[i - 1] /= h;
	    mean /= h;
	}
	else
	{
	    initMean[i - 1] = Eigen::VectorXd::Zero(Config::rank);
	    mean = Eigen::VectorXd::Zero(Utility::ndof);
	}
	std::ofstream output;
	output.open("mean" + std::to_string(counter) + "_" + std::to_string(i) + ".txt");
	output << mean.transpose() << std::endl;
	output.close();
    }

    std::vector<std::shared_ptr<const Sample>> minSamplesList = minSample->getMinSamplesList();
    std::vector<Eigen::VectorXd> minTrajectory;
    for (std::shared_ptr<const Sample> s: minSamplesList)
    {
	std::cout << s->cost << " ";
	minTrajectory.insert(minTrajectory.end(), s->trajectory.begin(), s->trajectory.end());
    }
    std::cout << std::endl;

    std::ofstream output;
    output.open(taskFileName + std::to_string(counter) + ".txt");
    for (const Eigen::VectorXd &v: Utility::bvhs[omp_get_thread_num()].frameToEulerAngle(minTrajectory))
	output << v.transpose() << std::endl;
    output.close();
    output.open("delta" + std::to_string(counter) + ".txt");
    for (std::shared_ptr<const Sample> s: minSamplesList)
	output << s->delta.transpose() << std::endl;
    output.close();
    size_t tRank = Config::rank;
    Config::rank += Config::dRank;
    if (Config::rank > Utility::ndof)
	Config::rank = Utility::ndof;
    Config::initSigma *= 0.7;
    for (VectorXd &m: initMean)
    {
	VectorXd v = m;
	m = VectorXd::Zero(Config::rank);
	m.head(tRank) = v;
    }
}
