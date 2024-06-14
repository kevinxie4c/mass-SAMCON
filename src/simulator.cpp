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
    if (Config::collisionDetector == "ode")
	solver->setCollisionDetector(detector = dart::collision::OdeCollisionDetector::create());
    else if (Config::collisionDetector == "fcl")
	solver->setCollisionDetector(detector = dart::collision::FCLCollisionDetector::create());
    else
	solver->setCollisionDetector(detector = dart::collision::DARTCollisionDetector::create());
    skeletonGroup = detector->createCollisionGroup(skeleton.get());
    floorGroup = detector->createCollisionGroup(floor.get());
}

void Simulator::setPose(const Eigen::VectorXd &pose, const Eigen::VectorXd &vel)
{
    skeleton->setPositions(pose);
    skeleton->setVelocities(vel);
}

#ifdef NDEBUG
bool Simulator::driveTo(const Eigen::VectorXd &ref, const std::vector<Eigen::VectorXd> &iforce, std::vector<Eigen::VectorXd> &resultTrajectory, std::vector<Eigen::VectorXd> &forces, bool useID)
#else
bool Simulator::driveTo(const Eigen::VectorXd &ref, const std::vector<Eigen::VectorXd> &iforce, std::vector<Eigen::VectorXd> &resultTrajectory, std::vector<Eigen::Vector3d> &com, std::vector<Eigen::Vector3d> &mmt, std::vector<Eigen::VectorXd> &forces, bool useID)
#endif
{
    for (size_t i = 0; i < Config::groupNum; ++i)
    {
	Eigen::Vector3d rightHipForce;
	if (Config::useCompensator)
	{
	    bool rightFootInAir = true;
	    dart::collision::CollisionResult result;
	    dart::collision::CollisionOption option;
	    if (skeletonGroup->collide(floorGroup.get(), option, &result))
	    {
		for (auto it: result.getCollidingBodyNodes())
		{
		    if (it->getName() == Utility::rightFootName)
		    {
			rightFootInAir = false;
			break;
		    }
		}
	    }
	    if (rightFootInAir)
	    {
		BodyNode *rightFoot = skeleton->getBodyNode(Utility::rightFootIndex);
		Eigen::MatrixXd jac = skeleton->getLinearJacobian(rightFoot, rightFoot->getLocalCOM());
		Eigen::MatrixXd J = jac.block(0, Utility::rightHipDofIdx, 3, 3);
		Eigen::Vector3d x = skeleton->getCOM() - rightFoot->getCOM();
		Eigen::Vector3d dx = rightFoot->getCOMLinearVelocity();
		x.y() = 0;
		dx.y() = 0;
		double l = sqrt(dx.norm() * dx.norm() * x.y() / Config::gravity);
		Eigen::Vector3d desire = l * dx.normalized();
		Eigen::Vector3d diff = desire + x; // desire - (-x);
		rightHipForce = J.transpose() * (Config::k_cmp * diff - Config::d_cmp * dx);
	    }
	}
	if (Config::stablePD)
	{
	    Eigen::VectorXd q = skeleton->getPositions();
	    Eigen::VectorXd dq = skeleton->getVelocities();
	    Eigen::MatrixXd invM = (skeleton->getMassMatrix() + Utility::mKd * skeleton->getTimeStep()).inverse();
	    Eigen::VectorXd p = -Utility::mKp * skeleton->getPositionDifferences(q + dq * skeleton->getTimeStep(), ref);
	    Eigen::VectorXd d = -Utility::mKd * dq;
	    Eigen::VectorXd qddot = invM * (-skeleton->getCoriolisAndGravityForces() + p + d + skeleton->getConstraintForces());
	    Eigen::VectorXd force = p + d -Utility::mKd * qddot * skeleton->getTimeStep();
	    if (useID)
	       	force += iforce[i];
	    force.segment(Utility::rightHipDofIdx, 3) += rightHipForce;
	    //Eigen::VectorXd force = iforce[i];
	    forces.push_back(force);
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
		force.segment(Utility::rightHipDofIdx, 3) += rightHipForce;
		forces.push_back(force);
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
