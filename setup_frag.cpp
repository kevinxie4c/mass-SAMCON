#include <cmath>
#include "mass.h"

using namespace std;
using namespace Eigen;

// set up frags
void setUpFrags(bool useMass)
{
    for (size_t i = 0; i < numFrag; ++i)
    {
	frags.push_back(ControlFragment());
	ControlFragment &f = frags.back();
	f.m = Eigen::MatrixXd::Zero(11, 18);
	f.sigma = Eigen::MatrixXd::Identity(11, 11) * Config::initSigma * Config::initSigma;
	f.a = Eigen::VectorXd::Zero(11);
	VectorXd pose;
	VectorXd vel;
	Utility::setStateAt(i * Config::groupNum + Config::groupNum, pose, vel);
	Utility::bvhs[omp_get_thread_num()].skeleton->setPositions(pose);
	Utility::bvhs[omp_get_thread_num()].skeleton->setVelocities(vel);
	std::vector<Eigen::Vector3d> q;
	std::vector<Eigen::Vector3d> omega;
	for (size_t i = 0; i < Utility::bvhs[omp_get_thread_num()].skeleton->getJoints().size(); ++i)
	{
	    // root has 6 DOFs, so?
	    Joint *joint = Utility::bvhs[omp_get_thread_num()].skeleton->getJoints()[i];
	    if (i == Utility::rootIndex)
	    {
		Eigen::VectorXd u = joint->getPositions();
		Eigen::Vector3d v(u[0], u[1], u[2]);
		q.push_back(v);
		u = joint->getVelocities();
		Eigen::Vector3d w(u[0], u[1], u[2]);
		omega.push_back(w);
	    }
	    else
	    {
		if (joint->getPositions().rows() == 3)
		{
		    q.push_back(joint->getPositions());
		    omega.push_back(joint->getVelocities());
		}
		else
		{
		    q.push_back(Eigen::Vector3d(joint->getPosition(0), 0, 0));
		    omega.push_back(Eigen::Vector3d(joint->getVelocity(0), 0, 0));
		}
	    }
	}
	f.endQ = q;
	f.endVel = omega;
	std::vector<Eigen::Vector3d> p;
	for (BodyNode const *bodyNode: Utility::bvhs[omp_get_thread_num()].skeleton->getBodyNodes())
	    p.push_back(bodyNode->getCOM());
	f.endPos = p;
	f.endCOM = Utility::bvhs[omp_get_thread_num()].skeleton->getCOM();
	f.endCOMVel = Utility::bvhs[omp_get_thread_num()].skeleton->getCOMLinearVelocity();

	f.tracked = pose;
	f.vel = vel;
	if (useMass)
	{
	    Eigen::MatrixXd mass = Utility::bvhs[omp_get_thread_num()].skeleton->getMassMatrix();
	    Eigen::EigenSolver<Eigen::MatrixXd> es(mass);
	    Eigen::MatrixXd B = es.eigenvectors().real();
	    Eigen::VectorXd D = es.eigenvalues().real();
	    for (int i = D.rows() - 2; i >= 0; --i)
		for (int j = 0; j <= i; ++j)
		    if (D[j] < D[j + 1])
		    {
			D.row(j).swap(D.row(j + 1));
			B.col(j).swap(B.col(j + 1));
		    }
	    // rank + (counter - 1) * dRank?
	    //f.transformation = B.leftCols(Config::rank + 6).rightCols(Config::rank) * Config::scaleMassMatrix;
	    double scaleMassMatrix = Config::scaleMassMatrix;
	    if (Config::autoScaleMassMatrix)
	    {
		//scaleMassMatrix = sqrt(Utility::ndof / mass.trace());
		scaleMassMatrix = sqrt(1.0 / D[6]);
		std::cout << "scaleMassMatrix " << i << ": " << scaleMassMatrix << std::endl;
	    }
	    if (Config::useEigenvalueScale)
		f.transformation = B * scaleMassMatrix * DiagonalMatrix<double, Dynamic, Dynamic>(D.array().sqrt().matrix());
	    else
		f.transformation = B * scaleMassMatrix;
	}
	else
	    f.transformation = Eigen::MatrixXd::Identity(Utility::ndof, Utility::ndof);
	
    }

    if (Config::useAFT)
    {
	Simulator sim(Utility::bvhs[omp_get_thread_num()]);
	// assume that we start from first frame
	VectorXd initPose, initVel;
	Utility::setStateAt(0, initPose, initVel);
	sim.setPose(initPose, initVel);
	vector<VectorXd> dummyTraj;
#ifndef NDEBUG
	vector<Vector3d> dummyCom, dummyMmt;
	vector<VectorXd> dummyForces;
	sim.driveTo(frags[0].tracked, dummyTraj, dummyCom, dummyMmt, dummyForces);
#else
	sim.driveTo(frags[0].tracked, dummyTraj);
#endif
	frags[0].aftOffset = Utility::bvhs[omp_get_thread_num()].skeleton->getPositionDifferences(frags[0].tracked, sim.skeleton->getPositions());

	for (size_t i = 1; i < frags.size(); ++i)
	{
	    sim.setPose(frags[i - 1].tracked, frags[i - 1].vel);
#ifndef NDEBUG
	    sim.driveTo(frags[i].tracked, dummyTraj, dummyCom, dummyMmt, dummyForces);
#else
	    sim.driveTo(frags[i].tracked, dummyTraj);
#endif
	    frags[i].aftOffset = Utility::bvhs[omp_get_thread_num()].skeleton->getPositionDifferences(frags[i].tracked, sim.skeleton->getPositions());
	}
    }
}
