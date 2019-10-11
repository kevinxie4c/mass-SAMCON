#include <algorithm>
#include <limits>
#include "simulator.h"
#include "utility.h"

Simulator::Simulator(BVHData &bvh)
{
    using namespace dart::simulation;
    world = World::create();
    world->setGravity(Eigen::Vector3d(0.0, Config::gravity, 0.0));
    world->setTimeStep(Config::timeStep);
    skeleton = bvh.skeleton->clone();
    world->addSkeleton(skeleton);
    world->addSkeleton(floor = Utility::createFloor());
    solver = world->getConstraintSolver();
    solver->setLCPSolver(dart::common::make_unique<dart::constraint::PGSLCPSolver>(world->getTimeStep()));
    solver->setCollisionDetector(detector = dart::collision::OdeCollisionDetector::create());
}

void Simulator::setPose(const Eigen::VectorXd &pose, const Eigen::VectorXd &vel)
{
    skeleton->setPositions(pose);
    skeleton->setVelocities(vel);
}

#ifdef NDEBUG
bool Simulator::driveTo(const Eigen::VectorXd &ref, const std::vector<Eigen::VectorXd> &iforce, std::vector<Eigen::VectorXd> &resultTrajectory)
#else
bool Simulator::driveTo(const Eigen::VectorXd &ref, const std::vector<Eigen::VectorXd> &iforce, std::vector<Eigen::VectorXd> &resultTrajectory, std::vector<Eigen::Vector3d> &com, std::vector<Eigen::Vector3d> &mmt, std::vector<Eigen::VectorXd> &forces)
#endif
{
    for (size_t i = 0; i < Config::groupNum; ++i)
    {
	if (Config::stablePD)
	{
	    Eigen::VectorXd q = skeleton->getPositions();
	    Eigen::VectorXd dq = skeleton->getVelocities();
	    Eigen::MatrixXd invM = (skeleton->getMassMatrix() + Utility::mKd * skeleton->getTimeStep()).inverse();
	    Eigen::VectorXd p = -Utility::mKp * skeleton->getPositionDifferences(q + dq * skeleton->getTimeStep(), ref);
	    Eigen::VectorXd d = -Utility::mKd * dq;
	    Eigen::VectorXd qddot = invM * (-skeleton->getCoriolisAndGravityForces() + p + d + skeleton->getConstraintForces());
	    Eigen::VectorXd force = p + d -Utility::mKd * qddot * skeleton->getTimeStep() + iforce[i];
	    //Eigen::VectorXd force = iforce[i];
#ifndef NDEBUG
	    forces.push_back(force);
#endif
	    for (size_t k = 0; k < Config::stepPerFrame; ++k)
	    {
		skeleton->setForces(force);
		world->step();
	    }
	}
	else
	{
	    for (size_t k = 0; k < Config::stepPerFrame; ++k)
	    {
		Eigen::VectorXd q = skeleton->getPositions();
		Eigen::VectorXd dq = skeleton->getVelocities();
		Eigen::VectorXd p = -Utility::mKp * skeleton->getPositionDifferences(q, ref);
		Eigen::VectorXd d = -Utility::mKd * dq;
		// forces.size() != groupNum. need to fix
		Eigen::VectorXd force = p + d + iforce[i];
#ifndef NDEBUG
		forces.push_back(force);
#endif
		skeleton->setForces(force);
		world->step();
	    }
	    
	}
	
	resultTrajectory.push_back(skeleton->getPositions());
	const BodyNode *bn = skeleton->getRootBodyNode();
#ifndef NDEBUG
	com.push_back(skeleton->getCOM());
	mmt.push_back(bn->getCOMLinearVelocity());
#endif
    }
    return true;
}
