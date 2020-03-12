#include "sample.h"

Sample::Sample(const Eigen::VectorXd &pose, const Eigen::VectorXd &vel): cost(0), parent(nullptr), resultPose(pose), resultVel(vel)
{
    if (Config::zeroInitialVelocities)
	resultVel = Eigen::VectorXd::Zero(resultPose.rows());
    trajectory.push_back(resultPose);
#ifndef NDEBUG
    com.push_back(Eigen::Vector3d::Zero());
    mmt.push_back(Eigen::Vector3d::Zero());
#endif
}

Sample::Sample(std::shared_ptr<Sample> parent, ControlFragment &cf, const Eigen::VectorXd &delta, const Eigen::VectorXd &kernel, Simulator &sim): parent(parent)
{
    sim.setPose(parent->resultPose, parent->resultVel);
    this->delta = delta;
    this->kernel = kernel;
    for (size_t i = 0; i < 6; ++i)
	this->delta[i] = 0;
    ref = cf.tracked + this->delta;
    if (Config::useAFT)
	ref = sim.skeleton->getPositionDifferences(ref, -cf.aftOffset);
#ifdef NDEBUG
    if (sim.driveTo(ref, cf.iforces, trajectory))
#else
    if (sim.driveTo(ref, cf.iforces, trajectory, com, mmt, forces))
#endif
    {
	resultPose = sim.skeleton->getPositions();
	resultVel = sim.skeleton->getVelocities();
	cost = Utility::costFunc(sim.skeleton, cf, zmp, et);
	totalCost = parent->totalCost + cost;
	accumCost = cost;
	if (parent != nullptr)
	    parent->update(height + 1, accumCost + parent->cost);
    }
    else
    {
	cost = NAN;
	std::cerr << "NAN" << std::endl;
    }
}

std::vector<Eigen::VectorXd> Sample::getTarget()
{
    std::vector<Eigen::VectorXd> list;
    //std::cout << cost << " ";
    list.push_back(ref);
    std::shared_ptr<const Sample> ptr = parent;
    //while (ptr != nullptr && ptr->parent != nullptr)
    while (ptr != nullptr)
    {
	//std::cout << ptr->cost << " ";
	list.push_back(ptr->ref);
	ptr = ptr->parent;
    }
    //std::cout << std::endl;
    std::reverse(std::begin(list), std::end(list));
    return list;
}

std::vector<std::shared_ptr<const Sample>> Sample::getMinSamplesList()
{
    std::vector<std::shared_ptr<const Sample>> list;
    //std::cout << cost << " ";
    list.push_back(shared_from_this());
    std::shared_ptr<const Sample> ptr = parent;
    //while (ptr != nullptr && ptr->parent != nullptr)
    while (ptr != nullptr)
    {
	//std::cout << ptr->cost << " ";
	list.push_back(ptr);
	ptr = ptr->parent;
    }
    //std::cout << std::endl;
    std::reverse(std::begin(list), std::end(list));
    return list;
}

std::vector<Eigen::VectorXd> Sample::getTrajectory() // should use getMinSamplesList instead
{
    std::vector<Eigen::VectorXd> list;
    //std::cout << cost << " ";
    for (auto it = trajectory.crbegin(); it != trajectory.crend(); ++it)
	list.push_back(*it);
    std::shared_ptr<const Sample> ptr = parent;
    //while (ptr != nullptr && ptr->parent != nullptr)
    while (ptr != nullptr)
    {
	//std::cout << ptr->cost << " ";
	for (auto it = ptr->trajectory.crbegin(); it != ptr->trajectory.crend(); ++it)
	    list.push_back(*it);
	ptr = ptr->parent;
    }
    //std::cout << std::endl;
    std::reverse(std::begin(list), std::end(list));
    return list;
}

std::vector<Eigen::VectorXd> Sample::getRef() // should use getMinSamplesList instead
{
    std::vector<Eigen::VectorXd> list;
    for (size_t i = 0; i < 2 * Config::groupNum; ++i)
	list.push_back(ref);
    std::shared_ptr<const Sample> ptr = parent;
    while (ptr != nullptr && ptr->parent != nullptr)
    {
	for (size_t i = 0; i < Config::groupNum; ++i)
	    list.push_back(ptr->ref);
	ptr = ptr->parent;
    }
    //std::cout << std::endl;
    std::reverse(std::begin(list), std::end(list));
    return list;
}

bool operator<(const Sample &lhs, const Sample &rhs)
{
    return lhs.cost < rhs.cost;
}

void Sample::update(size_t height, double accumCost)
{
    // can be optimized?
    if (height > this->height || (height == this->height && accumCost < this->accumCost))
    {
	this->height = height;
	this->accumCost = accumCost;
	if (parent != nullptr)
	    parent->update(height + 1, accumCost + parent->cost);
    }
}

bool CostCmp::operator()(const std::shared_ptr<Sample> &lhs, const std::shared_ptr<Sample> &rhs)
{
    return *lhs < *rhs;
}

bool SampleComparison::operator()(const std::shared_ptr<Sample> &lhs, const std::shared_ptr<Sample> &rhs)
{
    if (lhs->height != rhs->height)
	return lhs->height < rhs->height;
    else
	return rhs->accumCost < lhs->accumCost;
}
