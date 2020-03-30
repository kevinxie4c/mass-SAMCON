#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <dart/dart.hpp>
#include <dart/collision/ode/OdeCollisionDetector.hpp>
#include <Eigen/Core>
#include "config_loader.h"
#include "BVHData.h"
#include "utility.h"

class Simulator
{
    public:
	dart::simulation::WorldPtr world;
	dart::dynamics::SkeletonPtr skeleton;
	dart::dynamics::SkeletonPtr floor;
	dart::constraint::ConstraintSolver *solver;
	dart::collision::CollisionDetectorPtr detector;
	dart::collision::CollisionGroupPtr skeletonGroup, floorGroup;

	Simulator(BVHData &bvh);

	void setPose(const Eigen::VectorXd &pose, const Eigen::VectorXd &vel);

	const double largeNum = 100;
#ifdef NDEBUG
	bool driveTo(const Eigen::VectorXd &ref, const std::vector<Eigen::VectorXd> &iforces, std::vector<Eigen::VectorXd> &resultTrajectory);
#else
	bool driveTo(const Eigen::VectorXd &ref, const std::vector<Eigen::VectorXd> &iforces, std::vector<Eigen::VectorXd> &resultTrajectory, std::vector<Eigen::Vector3d> &com, std::vector<Eigen::Vector3d> &mmt, std::vector<Eigen::VectorXd> &forces);
#endif

};
#endif
