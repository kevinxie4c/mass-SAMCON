#ifndef SAMPLE_H
#define SAMPLE_H

#include <Eigen/Core>
#include "simulator.h"
#include "config_loader.h"
#include "utility.h"
#include "control_fragment.h"

class Sample: public std::enable_shared_from_this<Sample>
{
    public:
	double cost;
	Utility::ErrorTerms et;
	std::shared_ptr<Sample> parent = nullptr;
	Eigen::VectorXd delta;
	Eigen::VectorXd kernel;
	Eigen::VectorXd ref;
	Eigen::VectorXd resultPose;
	Eigen::VectorXd resultVel;
	Eigen::Vector3d zmp;
	double totalCost = 0; // total cost from ascendant to sample
	double accumCost = 0; // accumulative cost of the best path from the sample to the leaves
	size_t height = 0;
	std::vector<Eigen::VectorXd> trajectory;
#ifndef NDEBUG
	std::vector<Eigen::Vector3d> com;
	std::vector<Eigen::Vector3d> mmt;
	std::vector<Eigen::VectorXd> forces;
#endif

	//Sample(Eigen::VectorXd pose): cost(0), parent(nullptr), resultPose(pose), resultVel(Eigen::VectorXd::Zero(bvh.getChannelSize())) {}

	Sample(const Eigen::VectorXd &pose, const Eigen::VectorXd &vel);

	Sample(std::shared_ptr<Sample> parent, ControlFragment &cf, const Eigen::VectorXd &delta, const Eigen::VectorXd &kernel, Simulator &sim, bool useID);

	std::vector<Eigen::VectorXd> getTarget();

	std::vector<std::shared_ptr<const Sample>> getMinSamplesList();

	std::vector<Eigen::VectorXd> getTrajectory(); // should use getMinSamplesList instead

	std::vector<Eigen::VectorXd> getRef(); // should use getMinSamplesList instead

	friend bool operator<(const Sample &lhs, const Sample &rhs);

    private:
	void update(size_t height, double accumCost);
};

class CostCmp
{
    public:
	bool operator()(const std::shared_ptr<Sample> &lhs, const std::shared_ptr<Sample> &rhs);
};

class SampleComparison
{
    public:
	bool operator()(const std::shared_ptr<Sample> &lhs, const std::shared_ptr<Sample> &rhs);
};

#endif
