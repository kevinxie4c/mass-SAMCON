#include <dart/dart.hpp>
#include <dart/collision/ode/OdeCollisionDetector.hpp>
#include <vector>
#include <queue>
#include <string>
#include <algorithm>
#include <cmath>
#include <memory>
#include <cassert>
#include <cfloat>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include "eigenmvn.h"
#include "BVHData.h"
#include "cmaes.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#endif

using namespace dart::dynamics;
using namespace dart::simulation;
using namespace dart::collision;

double default_stiffness = 500.0;
double default_damping = 50.0;
double frictionCoeff = 0.5;
double restitutionCoeff = 0.0;
double gravity = -9.8;
double ERP = 0.2;
double CFM = 0.01;
Eigen::MatrixXd mKp, mKd;

double wp= 5, wr = 3, we = 30, wb = 10;

BVHData bvh;
std::string bvhFileName = "walk.bvh";
std::string outputFileName = "result.txt";
size_t startFrame = 32;
size_t endFrame = 300;
double groundOffset = -0.03;
size_t sampleNum = 1200, saveNum = 400;
size_t groupNum = 5;
size_t slidingWindow = 50;
size_t trialMin = 5, trialMax = 20;
size_t notImproveMax = 5;
double goodEnough = 5.0, failThreshold = 20.0;
double frameTime = 0.008333;
double timeStep = 0.001;
//double timeStep = 0.008333;
double init_cov = 0.008;
size_t rank = 10;
double scaleMassMatrix = 0.3;
std::string cov_fileName = "default_cov.txt";
std::string stiffness_fileName = "default_stiffness.txt";
std::string damping_fileName = "default_damping.txt";
std::string mass_fileName = "default_mass.txt";
bool stablePD = true;
bool useSampleNumAsLambda = false;
double init_sigma = 0.1;
std::vector<double> init_cov_list;
std::vector<double> stiffness_list;
std::vector<double> damping_list;
std::vector<double> mass_list;
size_t stepPerFrame = frameTime / timeStep;
size_t rootIndex;
std::vector<size_t> endEffectorIndex;
std::string rootName = "Hip_joint";
std::vector<std::string> endEffectorName = { "Wrist", "Head", "Ankle" };

SkeletonPtr createFloor()
{
    SkeletonPtr floor = Skeleton::create("floor");

    // Give the floor a body
    BodyNodePtr body = floor->createJointAndBodyNodePair<WeldJoint>(nullptr).second;
    body->setFrictionCoeff(frictionCoeff);
    body->setRestitutionCoeff(restitutionCoeff);

    // Give the body a shape
    double floor_width = 20.0;
    double floor_height = 0.01;
    std::shared_ptr<BoxShape> box(
	    new BoxShape(Eigen::Vector3d(floor_width, floor_height, floor_width)));
    auto shapeNode = body->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(box);
    shapeNode->getVisualAspect()->setColor(dart::Color::Black());

    // Put the body into position
    Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
    tf.translation() = Eigen::Vector3d(0.0, groundOffset - floor_height / 2, 0.0);
    body->getParentJoint()->setTransformFromParentBodyNode(tf);

    return floor;
}

SkeletonPtr createBox()
{
    SkeletonPtr object = Skeleton::create("box");

    // Give the floor a body
    BodyNodePtr body = object->createJointAndBodyNodePair<FreeJoint>(nullptr).second;
    body->setFrictionCoeff(frictionCoeff);
    body->setRestitutionCoeff(restitutionCoeff);

    // Give the body a shape
    double box_width = 0.2;
    std::shared_ptr<BoxShape> box(
	    new BoxShape(Eigen::Vector3d(box_width, box_width, box_width)));
    auto shapeNode = body->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(box);
    shapeNode->getVisualAspect()->setColor(dart::Color::Black());

    // Put the body into position
    Eigen::Vector6d positions(Eigen::Vector6d::Zero());
    positions[4] = 0.5;
    object->getJoint(0)->setPositions(positions);
    return object;
}
std::vector<std::vector<Eigen::Vector3d>> target_q;
std::vector<std::vector<Eigen::Vector3d>> target_omega;
std::vector<std::vector<Eigen::Vector3d>> target_p;
std::vector<Eigen::Vector3d> target_COM;
std::vector<Eigen::Vector3d> target_COM_vel;

// set position and velocity at index
void setStateAt(size_t index, Eigen::VectorXd &pose, Eigen::VectorXd &vel)
{
    static SkeletonPtr skeleton = bvh.skeleton->clone();
    static SkeletonPtr t_skeleton = bvh.skeleton->clone(); // for velocity computation
    skeleton->setPositions(bvh.frame[index]);
    t_skeleton->setPositions(bvh.frame[index + 1]);
    // also need to set velocities
    // what's the rule when we use rotation vector?
    // assume that the first one is the root joint
    {
	// root joint
	Joint *joint = skeleton->getJoint(0);
	Eigen::Vector6d q = joint->getPositions();
	Eigen::Vector3d p = q.tail<3>(); // assume that last 3 is translation
	Eigen::Vector3d omega = q.head<3>(); // assume that first 3 is rotation, and rotation is in rotation vector form
	Joint *t_joint = t_skeleton->getJoint(0);
	Eigen::Vector6d t_q = t_joint->getPositions();
	Eigen::Vector3d t_p = t_q.tail<3>();
	Eigen::Vector3d t_omega = t_q.head<3>();
	Eigen::AngleAxisd r(omega.norm(), omega.normalized());
	Eigen::AngleAxisd t_r(t_omega.norm(), t_omega.normalized());
	Eigen::AngleAxisd diff_r(t_r * r.inverse());
	Eigen::Vector6d v;
	// rotation, translation
	//v << diff_r.axis() * (diff_r.angle() / frameTime), (t_p - p) / frameTime;
	v << 0, 0, 0, (t_p - p) / frameTime;
	joint->setVelocities(v);
    }
    // other joints
    for (size_t i = 1; i < skeleton->getJoints().size(); ++i)
    {
	Joint *joint = skeleton->getJoint(i);
	Eigen::Vector3d omega = joint->getPositions();
	Joint *t_joint = t_skeleton->getJoint(i);
	Eigen::Vector3d t_omega = t_joint->getPositions();
	Eigen::AngleAxisd r(omega.norm(), omega.normalized());
	Eigen::AngleAxisd t_r(t_omega.norm(), t_omega.normalized());
	Eigen::AngleAxisd diff_r(t_r * r.inverse());
	joint->setVelocities(diff_r.axis() * (diff_r.angle() / frameTime));
    }
    pose = skeleton->getPositions();
    vel = skeleton->getVelocities();
}

void setMassFor(BVHData bvh)
{
    size_t i = 0;
    for (BodyNode *bn: bvh.skeleton->getBodyNodes())
	bn->setMass(mass_list[i++]);
}

void initCost(std::string filename)
{
    BVHData target;
    target.loadBVH(filename);
    setMassFor(target);
    for (size_t index = 0; index < target.frame.size(); ++index)
    {
	Eigen::VectorXd pose;
	Eigen::VectorXd vel;
	if (index < target.frame.size() - 1)
	    setStateAt(index, pose, vel);
	else
	{
	    // use velocity at index - 1 for last frame
	    setStateAt(index - 1, pose, vel);
	    pose = target.frame.back();
	}
	target.skeleton->setPositions(pose);
	target.skeleton->setVelocities(vel);
	std::vector<Eigen::Vector3d> q;
	std::vector<Eigen::Vector3d> omega;
	for (size_t i = 0; i < target.skeleton->getJoints().size(); ++i)
	{
	    // root has 6 DOFs, so?
	    Joint *joint = target.skeleton->getJoints()[i];
	    if (i == rootIndex)
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
		q.push_back(joint->getPositions());
		omega.push_back(joint->getVelocities());
	    }
	}
	target_q.push_back(q);
	target_omega.push_back(omega);
	std::vector<Eigen::Vector3d> p;
	for (BodyNode const *bodyNode: target.skeleton->getBodyNodes())
	    p.push_back(bodyNode->getCOM());
	target_p.push_back(p);
	target_COM.push_back(target.skeleton->getCOM());
	target_COM_vel.push_back(target.skeleton->getCOMLinearVelocity());
    }
}

double costFunc(const SkeletonPtr skeleton, size_t index)
{
    static size_t n = skeleton->getJoints().size();
    double err_p = 0;
    std::vector<Eigen::Vector3d> qlist = target_q[index];
    std::vector<Eigen::Vector3d> omegalist = target_omega[index];
    for (size_t i = 0; i < skeleton->getJoints().size(); ++i)
    {
	const Joint *joint = skeleton->getJoints()[i];
	Eigen::Vector3d q;
	Eigen::Vector3d qr = qlist[i];
	Eigen::Vector3d omega;
	Eigen::Vector3d omega_r = omegalist[i];
	if (i == rootIndex)
	{
	    Eigen::VectorXd u = joint->getPositions();
	    q = Eigen::Vector3d(u[0], u[1], u[2]);
	    u = joint->getVelocities();
	    omega = Eigen::Vector3d(u[0], u[1], u[2]);
	}
	else
	{
	    q = joint->getPositions();
	    omega= joint->getVelocities();
	}
	Eigen::Quaterniond quat1(Eigen::AngleAxisd(q.norm(), q.normalized()));
	Eigen::Quaterniond quat2(Eigen::AngleAxisd(qr.norm(), qr.normalized()));
	err_p += fabs(quat1.angularDistance(quat2)) + (omega - omega_r).norm() * 0.1;
    }
    err_p /= n;

    double err_r = 0;
    {
	const Joint *joint = skeleton->getJoints()[rootIndex];
	Eigen::VectorXd u = joint->getPositions();
	Eigen::Vector3d q(u[0], u[1], u[2]);
	Eigen::Vector3d qr = qlist[rootIndex];
	Eigen::Quaterniond quat1(Eigen::AngleAxisd(q.norm(), q.normalized()));
	Eigen::Quaterniond quat2(Eigen::AngleAxisd(qr.norm(), qr.normalized()));
	err_r += fabs(quat1.angularDistance(quat2));
	u = joint->getVelocities();
	Eigen::Vector3d omega(u[0], u[1], u[2]);
	Eigen::Vector3d omega_r = omegalist[rootIndex];
	err_r += (omega - omega_r).norm() * 0.1;
    }

    static size_t k = endEffectorName.size();
    double err_e = 0;
    std::vector<Eigen::Vector3d> plist = target_p[index];
    for (size_t i: endEffectorIndex)
    {
	const BodyNode *node = skeleton->getBodyNodes()[i];
	Eigen::Vector3d p = node->getCOM();
	Eigen::Vector3d pr = plist[i];
	err_e += fabs(p.y() - pr.y());
	//err_e += (p - pr).norm();
    }
    err_e /= k;

    double err_b = 0;
    Eigen::Vector3d COM = skeleton->getCOM();
    Eigen::Vector3d COMr = target_COM[index];
    for (size_t i: endEffectorIndex)
    {
	const BodyNode *node = skeleton->getBodyNodes()[i];
	Eigen::Vector3d p = node->getCOM();
	Eigen::Vector3d pr = plist[i];
	Eigen::Vector3d rci = COM - p;
	rci.y() = 0;
	Eigen::Vector3d rci_r = COMr - pr;
	rci_r.y() = 0;
	err_b += (rci - rci_r).norm();
    }
    static double h = 1.6;
    err_b /= h * k;
    err_b += (skeleton->getCOMLinearVelocity() - target_COM_vel[index]).norm() * 0.1;

    return wp * err_p + wr * err_r + we * err_e + wb * err_b;
}

class Simulator
{
    public:
	WorldPtr world;
	SkeletonPtr skeleton;
	Simulator()
	{
	    world = World::create();
	    dart::constraint::ConstraintSolver *solver = world->getConstraintSolver();
	    solver->setLCPSolver(dart::common::make_unique<dart::constraint::PGSLCPSolver>(world->getTimeStep()));
	    solver->setCollisionDetector(dart::collision::OdeCollisionDetector::create());
	    world->setGravity(Eigen::Vector3d(0.0, gravity, 0.0));
	    world->setTimeStep(timeStep);
	    skeleton = bvh.skeleton->clone();
	    world->addSkeleton(skeleton);
	    world->addSkeleton(createFloor());
	}

	void setPose(Eigen::VectorXd pose, Eigen::VectorXd vel)
	{
	    skeleton->setPositions(pose);
	    skeleton->setVelocities(vel);
	}

	const double largeNum = 100;
	bool driveTo(Eigen::VectorXd ref, std::vector<Eigen::VectorXd> &resultTrajectory)
	{
	    for (size_t i = 0; i < groupNum; ++i)
	    {
		if (stablePD)
		{
		    Eigen::VectorXd q = skeleton->getPositions();
		    Eigen::VectorXd dq = skeleton->getVelocities();
		    Eigen::MatrixXd invM = (skeleton->getMassMatrix() + mKd * skeleton->getTimeStep()).inverse();
		    Eigen::VectorXd p = -mKp * (q + dq * skeleton->getTimeStep() - ref);
		    Eigen::VectorXd d = -mKd * dq;
		    Eigen::VectorXd qddot = invM * (-skeleton->getCoriolisAndGravityForces() + p + d + skeleton->getConstraintForces());
		    Eigen::VectorXd force = p + d - mKd * qddot * skeleton->getTimeStep();
		    for (size_t k = 0; k < stepPerFrame; ++k)
		    {
			skeleton->setForces(force);
			world->step();
		    }
		}
		else
		{
		    for (size_t k = 0; k < stepPerFrame; ++k)
		    {
			Eigen::VectorXd q = skeleton->getPositions();
			Eigen::VectorXd dq = skeleton->getVelocities();
			Eigen::VectorXd p = -mKp * (q - ref);
			Eigen::VectorXd d = -mKd * dq;
			Eigen::VectorXd force = p + d;
			skeleton->setForces(force);
			world->step();
		    }
		}
		resultTrajectory.push_back(skeleton->getPositions());
	    }
	    return true;
	}
};

std::vector<Simulator> simulators;

class Sample
{
    public:
	double cost;
	Eigen::VectorXd delta;
	Eigen::VectorXd kernel;
	Eigen::VectorXd ref;
	Eigen::VectorXd resultPose;
	Eigen::VectorXd resultVel;
	std::shared_ptr<Sample> parent = nullptr;
	double totalCost = 0; // total cost from ascendant to sample
	double accumCost = 0; // accumulative cost of the best path from the sample to the leaves
	size_t height = 0;
	std::vector<Eigen::VectorXd> trajectory;

	//Sample(Eigen::VectorXd pose): cost(0), parent(nullptr), resultPose(pose), resultVel(Eigen::VectorXd::Zero(bvh.getChannelSize())) {}

	Sample(size_t index): cost(0), parent(nullptr)
	{
	    setStateAt(index, resultPose, resultVel);
	    trajectory.push_back(resultPose);
	}

	Sample(std::shared_ptr<Sample> parent, size_t index, Eigen::VectorXd delta, Eigen::VectorXd kernel, size_t simIndex): parent(parent)
	{
	    Simulator &sim = simulators[simIndex];
	    sim.setPose(parent->resultPose, parent->resultVel);
	    this->delta = delta;
	    this->kernel = kernel;
	    for (size_t i = 0; i < 6; ++i)
		delta[i] = 0;
	    ref = bvh.frame[index] + delta;
	    //std::cout << sim.skeleton->getJoint(0)->getPositions() << std::endl;
	    if (sim.driveTo(ref, trajectory))
	    {
		resultPose = sim.skeleton->getPositions();
		resultVel = sim.skeleton->getVelocities();
		cost = costFunc(sim.skeleton, index);
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

	std::vector<Eigen::VectorXd> getTarget()
	{
	    std::vector<Eigen::VectorXd> list;
	    std::cout << cost << " ";
	    list.push_back(ref);
	    std::shared_ptr<const Sample> ptr = parent;
	    while (ptr != nullptr && ptr->parent != nullptr)
	    {
		std::cout << ptr->cost << " ";
		list.push_back(ptr->ref);
		ptr = ptr->parent;
	    }
	    std::cout << std::endl;
	    std::reverse(std::begin(list), std::end(list));
	    return list;
	}

	std::vector<Eigen::VectorXd> getTrajectory()
	{
	    std::vector<Eigen::VectorXd> list;
	    std::cout << cost << " ";
	    for (auto it = trajectory.crbegin(); it != trajectory.crend(); ++it)
		list.push_back(*it);
	    std::shared_ptr<const Sample> ptr = parent;
	    while (ptr != nullptr && ptr->parent != nullptr)
	    {
		std::cout << ptr->cost << " ";
		for (auto it = ptr->trajectory.crbegin(); it != ptr->trajectory.crend(); ++it)
		    list.push_back(*it);
		ptr = ptr->parent;
	    }
	    std::cout << std::endl;
	    std::reverse(std::begin(list), std::end(list));
	    return list;
	}

	friend bool operator<(const Sample &lhs, const Sample &rhs)
	{
	    return lhs.cost < rhs.cost;
	}

    private:
	void update(size_t height, double accumCost)
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
};

class CostCmp
{
    public:
	bool operator()(const std::shared_ptr<Sample> &lhs, const std::shared_ptr<Sample> &rhs)
	{
	    return *lhs < *rhs;
	}
};

class SampleComparison
{
    public:
	bool operator()(const std::shared_ptr<Sample> &lhs, const std::shared_ptr<Sample> &rhs)
	{
	    if (lhs->height != rhs->height)
		return lhs->height < rhs->height;
	    else
		return rhs->accumCost < lhs->accumCost;
	}
};

void simulate(const BVHData &bvh, std::vector<Eigen::VectorXd> reference, std::vector<Eigen::VectorXd> &result, std::vector<std::vector<Eigen::Vector3d>> &force)
{
    WorldPtr world = World::create();
    dart::constraint::ConstraintSolver *solver = world->getConstraintSolver();
    solver->setLCPSolver(dart::common::make_unique<dart::constraint::PGSLCPSolver>(world->getTimeStep()));
    solver->setCollisionDetector(dart::collision::OdeCollisionDetector::create());
    world->setGravity(Eigen::Vector3d(0.0, gravity, 0.0));
    world->setTimeStep(timeStep);
    SkeletonPtr skeleton = bvh.skeleton->clone();
    world->addSkeleton(skeleton);
    world->addSkeleton(createFloor());

    // also need to set velocities
    // what's the rule when we use rotation vector?
    // assume that the first one is the root joint
    Eigen::VectorXd initPose;
    Eigen::VectorXd initVel;
    setStateAt(startFrame, initPose, initVel);
    skeleton->setPositions(initPose);
    skeleton->setVelocities(initVel);

    result.push_back(bvh.frame[startFrame]);
    for (size_t i = 0; i < reference.size(); ++i)
    {
	Eigen::VectorXd ref = reference[i];
	for (size_t j = 0; j < groupNum; ++j)
	{
	    if (stablePD)
	    {
		Eigen::VectorXd q = skeleton->getPositions();
		Eigen::VectorXd dq = skeleton->getVelocities();
		Eigen::MatrixXd invM = (skeleton->getMassMatrix() + mKd * skeleton->getTimeStep()).inverse();
		Eigen::VectorXd p = -mKp * (q + dq * skeleton->getTimeStep() - ref);
		Eigen::VectorXd d = -mKd * dq;
		Eigen::VectorXd qddot = invM * (-skeleton->getCoriolisAndGravityForces() + p + d + skeleton->getConstraintForces());
		Eigen::VectorXd force = p + d - mKd * qddot * skeleton->getTimeStep();
		for (size_t k = 0; k < stepPerFrame; ++k)
		{
		    skeleton->setForces(force);
		    world->step();
		}
	    }
	    else
	    {
		for (size_t k = 0; k < stepPerFrame; ++k)
		{
		    Eigen::VectorXd q = skeleton->getPositions();
		    Eigen::VectorXd dq = skeleton->getVelocities();
		    Eigen::VectorXd p = -mKp * (q - ref);
		    Eigen::VectorXd d = -mKd * dq;
		    Eigen::VectorXd force = p + d;
		    skeleton->setForces(force);
		    world->step();
		}
	    }
	    result.push_back(skeleton->getPositions());
	    CollisionResult cr = world->getLastCollisionResult();
	    std::vector<Eigen::Vector3d> f;
	    for (const Contact c: cr.getContacts())
	    {
		f.push_back(c.point);
		f.push_back(c.force);
	    }
	    force.push_back(f);
	}
	std::cout << costFunc(skeleton, startFrame + (i + 1) * groupNum) << std::endl;
    }
}

std::vector<Eigen::VectorXd> getTarget(const BVHData &bvh)
{
    static size_t counter = 0;
    static std::vector<Eigen::VectorXd> init_mean;
    static std::vector<Eigen::MatrixXd> transformation;
    ++counter;
    std::ofstream f_cov;
    std::ofstream f_mean;
    f_cov.open(std::string("cov_") + std::to_string(counter) + ".txt");
    f_mean.open(std::string("mean_") + std::to_string(counter) + ".txt");
    size_t actualEnd = std::min(bvh.frame.size(), endFrame);
    size_t dim = bvh.getChannelSize();
    size_t vectorSize = (actualEnd - startFrame) / groupNum + 1;
    if (counter == 1)
    {
	for (size_t i = 0; i < vectorSize; ++i)
	    init_mean.push_back(Eigen::VectorXd::Zero(rank));
    }
    size_t trial = 0;
    std::vector<WeirdCMAES> cmaes;
    for (size_t i = 0; i < vectorSize; ++i)
    {
	bvh.setPositionAt(startFrame + i * groupNum);
	Eigen::MatrixXd mass = bvh.skeleton->getMassMatrix();
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
	if (counter == 1)
	    transformation.push_back(B.leftCols(rank) * D.head(rank).array().sqrt().matrix().asDiagonal() * scaleMassMatrix);
	cmaes.push_back(WeirdCMAES(rank, sampleNum, saveNum, init_sigma, init_mean[i]));
    }
    std::vector<size_t> generation(vectorSize, 0);
    std::vector<size_t> notImprove(vectorSize, 0);
    std::vector<double> minCost(vectorSize, DBL_MAX);
    std::cout << "_OPENMP " << _OPENMP << std::endl;
    std::cout << "omp_get_max_threads() = " << omp_get_max_threads() << std::endl;
    for (size_t i = 0; i < omp_get_max_threads(); ++i)
	simulators.push_back(Simulator());
    std::vector<std::vector<std::shared_ptr<Sample>>> savedSamples;
    std::vector<std::vector<std::shared_ptr<Sample>>> backupSamples;
    double backupMin = DBL_MAX;
    std::vector<std::shared_ptr<Sample>> tmpList;
    // also need to set velocities
    // what's the rule when we use rotation vector?
    // assume that the first one is the root joint
    for (size_t i = 0; i < saveNum; ++i)
	tmpList.push_back(std::make_shared<Sample>(startFrame));
    //tmpList.push_back(Sample(bvh.frame[startFrame]));
    savedSamples.push_back(tmpList);
    size_t i_begin = 1, i_end;
    while (startFrame + i_begin * (groupNum + 1) < actualEnd)
    {
	i_end = i_begin + slidingWindow;
	if (startFrame + i_end * groupNum > actualEnd)
	    i_end = (actualEnd - startFrame) / groupNum;
	for (size_t i = i_begin; i < i_end && startFrame + i * groupNum < actualEnd; ++i)
	{
	    std::cout << startFrame + i * groupNum << std::endl;
	    std::vector<std::shared_ptr<Sample>> &samples = savedSamples[i - 1];
	    std::priority_queue<std::shared_ptr<Sample>, std::vector<std::shared_ptr<Sample>>, CostCmp> queue;
#pragma omp parallel for	
	    for (size_t j = 0; j < samples.size(); ++j)
	    {
		std::shared_ptr<Sample> &sample = samples[j];
		for (size_t k = 0; k < sampleNum / samples.size(); ++k)
		{
		    Eigen::VectorXd kernel = cmaes[i].getSample();
		    Eigen::VectorXd delta = transformation[i] * kernel;
		    std::shared_ptr<Sample> t = std::make_shared<Sample>(sample, startFrame + i * groupNum, delta, kernel, omp_get_thread_num());
#pragma omp critical (queue_section)
		    {
			queue.push(t);
			if (queue.size() > saveNum)
			    queue.pop();
		    }
		}
	    }
	    std::vector<std::shared_ptr<Sample>> tmp;
	    while (!queue.empty())
	    {
		const std::shared_ptr<Sample> &t = queue.top();
		if (!std::isnan(t->cost))
		    tmp.push_back(queue.top());
		queue.pop();
	    }
	    std::cout << tmp.back()->cost << std::endl;
	    if (tmp.back()->cost > failThreshold)
	    {
		i_end = i;
		std::cout << "fail" << std::endl;
		break;
	    }
	    if (tmp.back()->cost < minCost[i])
	    {
		minCost[i] = tmp.back()->cost;
		notImprove[i] = 0;
	    }
	    else
		++notImprove[i];
	    if (i < savedSamples.size())
		savedSamples[i] = tmp;
	    else
		savedSamples.push_back(tmp);
	}
	std::shared_ptr<Sample> minSample = nullptr;
	double min = DBL_MAX;
	for (std::shared_ptr<Sample> &sample: savedSamples[i_end - 1])
	{
	    if (sample->totalCost < min)
	    {
		min = sample->totalCost;
		minSample = sample;
	    }
	}
	std::ofstream output;
	output.open(outputFileName + std::to_string(counter) + "_" + std::to_string(++trial) + ".txt");
	std::cout << "trial: " << counter << " - " << trial << std::endl;
	for (const Eigen::VectorXd &v: bvh.frameToEulerAngle(minSample->getTrajectory()))
	    output << v.transpose() << std::endl;
	output.close();
	if (i_end > backupSamples.size() || (i_end == backupSamples.size() && min < backupMin))
	{
	    for (size_t i = i_begin; i < i_end; ++i)
	    {
		std::vector<std::shared_ptr<Sample>> &tmp = savedSamples[i];
		std::priority_queue<std::shared_ptr<Sample>, std::vector<std::shared_ptr<Sample>>, SampleComparison> tQueue(tmp.cbegin(), tmp.cend());
		Eigen::MatrixXd X(rank, saveNum);
		for (size_t j = 0; j < saveNum; ++j)
		{
		    X.col(j) = tQueue.top()->kernel;
		    tQueue.pop();
		}
		cmaes[i].update(X);
	    }
	    backupSamples.clear();
	    for (size_t i = 0; i < i_end; ++i)
		backupSamples.push_back(savedSamples[i]);
	    backupMin = min;
	}
	else
	{
	    savedSamples.clear();
	    for (auto i: backupSamples)
		savedSamples.push_back(i);
	}
	std::cout << std::endl;
	for (size_t i = 0; i < vectorSize; ++i)
	{
	    generation[i] = cmaes[i].counteval / cmaes[i].lambda;
	    std::cout << generation[i] << " ";
	}
	std::cout << std::endl;
	for (size_t i = 0; i < vectorSize; ++i)
	    std::cout << cmaes[i].sigma << " ";
	std::cout << std::endl;
	for (size_t i = 0; i < vectorSize; ++i)
	{
	    f_cov << cmaes[i].D.prod() * cmaes[i].sigma << " ";
	    f_mean << cmaes[i].xmean.transpose() << std::endl;
	}
	f_cov << std::endl;
	f_mean << std::endl;
	/*
	   std::vector<Eigen::VectorXd> frame;
	   std::vector<std::vector<Eigen::Vector3d>> force;
	   std::vector<Eigen::VectorXd> target;
	   target.push_back(bvh.frame[startFrame]);
	   std::vector<Eigen::VectorXd> tmp = minSample->getTarget();
	   for (const Eigen::VectorXd &v: tmp)
	   {
	   for (size_t i = 0; i < groupNum; ++i)
	   target.push_back(v);
	   }
	   simulate(bvh, tmp, frame, force);
	   for (const Eigen::VectorXd &v: bvh.frameToEulerAngle(frame))
	   output << v.transpose() << std::endl;
	   output.close();
	   output.open(outputFileName + std::to_string(trial) + ".force.txt");
	   output << force.size() << std::endl;
	   for (const std::vector<Eigen::Vector3d> &list: force)
	   {
	   output << list.size() << std::endl;
	   for (const Eigen::Vector3d &v: list)
	   output << v << std::endl;
	   }
	   output.close();
	   output.open(outputFileName + std::to_string(trial) + ".target.txt");
	   for (const Eigen::VectorXd &v: bvh.frameToEulerAngle(target))
	   output << v.transpose() << std::endl;
	   output.close();
	 */
	if (generation[i_begin] < trialMin) continue;
	while ((startFrame + i_begin * (groupNum + 1) < actualEnd) && (generation[i_begin] > trialMax || notImprove[i_begin] > notImproveMax || minCost[i_begin] < goodEnough))
	{
	    ++i_begin;
	}
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
	init_mean[i] = Eigen::VectorXd::Zero(rank);
	for (const std::shared_ptr<Sample> &s: savedSamples[i])
	{
	    if (s->height >= 4)
	    {
		init_mean[i] += s->kernel * s->height;
		h += s->height;
	    }
	}
	if (h > 0)
	    init_mean[i] /= h;
	else
	    init_mean[i] = Eigen::VectorXd::Zero(rank);
    }
    init_sigma *= 0.7;
    f_cov.close();
    f_mean.close();

    return minSample->getTarget();
}

void setParameter(std::string parameter, std::string value)
{
    if (parameter == "bvhFileName")
	bvhFileName = value;
    if (parameter == "outputFileName")
	outputFileName = value;
    else if (parameter == "startFrame")
	startFrame = std::stoi(value);
    else if (parameter == "endFrame")
	endFrame = std::stoi(value);
    else if (parameter == "groundOffset")
	groundOffset = std::stod(value);
    else if (parameter == "sampleNum")
	sampleNum = std::stoi(value);
    else if (parameter == "saveNum")
	saveNum = std::stoi(value);
    else if (parameter == "groupNum")
	groupNum = std::stoi(value);
    else if (parameter == "frameTime")
	frameTime = std::stod(value);
    else if (parameter == "timeStep")
	timeStep = std::stod(value);
    else if (parameter == "init_cov")
	init_cov = std::stod(value);
    else if (parameter == "default_stiffness")
	default_stiffness = std::stod(value);
    else if (parameter == "default_damping")
	default_damping = std::stod(value);
    else if (parameter == "frictionCoeff")
	frictionCoeff = std::stod(value);
    else if (parameter == "restitutionCoeff")
	restitutionCoeff = std::stod(value);
    else if (parameter == "wp")
	wp = std::stod(value);
    else if (parameter == "wr")
	wr = std::stod(value);
    else if (parameter == "we")
	we = std::stod(value);
    else if (parameter == "wb")
	wb = std::stod(value);
    else if (parameter == "cov_fileName")
	cov_fileName = value;
    else if (parameter == "stiffness_fileName")
	stiffness_fileName = value;
    else if (parameter == "damping_fileName")
	damping_fileName = value;
    else if (parameter == "mass_fileName")
	mass_fileName = value;
    else if (parameter == "gravity")
	gravity = std::stod(value);
    else if (parameter == "stablePD")
	stablePD = std::stoi(value);
    else if (parameter == "ERP")
	ERP = std::stod(value);
    else if (parameter == "CFM")
	CFM = std::stod(value);
    else if (parameter == "slidingWindow")
	slidingWindow = std::stoi(value);
    else if (parameter == "useSampleNumAsLambda")
	useSampleNumAsLambda = std::stoi(value);
    else if (parameter == "init_sigma")
	init_sigma = std::stod(value);
    else if (parameter == "failThreshold")
	failThreshold = std::stod(value);
    else if (parameter == "rank")
	rank = std::stoi(value);
    else if (parameter == "scaleMassMatrix")
	scaleMassMatrix = std::stod(value);
}

std::vector<double> readVectorDoubleFrom(std::string fileName)
{
    std::ifstream input;
    input.open(fileName);
    std::vector<double> result;
    double d;
    while (input >> d)
	result.push_back(d);
    return result;
    input.close();
}


int main(int argc, char* argv[])
{
    if (argc != 1)
    {
	std::ifstream taskFile;
	taskFile.open(argv[1]);
	std::string parameter;
	std::string value;
	while (taskFile >> parameter >> value)
	    setParameter(parameter, value);
	taskFile.close();
    }
    dart::constraint::ContactConstraint::setErrorReductionParameter(ERP);
    dart::constraint::ContactConstraint::setConstraintForceMixing(CFM);

    init_cov_list = readVectorDoubleFrom(cov_fileName);
    stiffness_list = readVectorDoubleFrom(stiffness_fileName);
    damping_list = readVectorDoubleFrom(damping_fileName);
    mass_list = readVectorDoubleFrom(mass_fileName);
    stepPerFrame = frameTime / timeStep;

    bvh.loadBVH(bvhFileName);
    setMassFor(bvh);

    for (size_t i = 0; i < bvh.skeleton->getJoints().size(); ++i)
    {
	std::cout << bvh.skeleton->getJoints()[i]->getName() << std::endl;
	if (bvh.skeleton->getJoints()[i]->getName() == rootName)
	{
	    rootIndex = i;
	    break;
	}
    }
    std::cout << rootIndex << std::endl;
    for (size_t i = 0; i < bvh.skeleton->getBodyNodes().size(); ++i)
    {
	std::cout << bvh.skeleton->getBodyNodes()[i]->getName() << std::endl;
	for (const std::string& str: endEffectorName)
	    if (bvh.skeleton->getBodyNodes()[i]->getName().find(str) != std::string::npos)
	    {
		endEffectorIndex.push_back(i);
		break;
	    }
    }
    for (int i: endEffectorIndex)
	std::cout << i << std::endl;
    initCost(bvhFileName);

    size_t numDofs = bvh.skeleton->getDofs().size();
    mKp = Eigen::MatrixXd::Zero(numDofs, numDofs);
    mKd = Eigen::MatrixXd::Zero(numDofs, numDofs);

    for (size_t i = 0; i < numDofs; ++i)
    {
	mKp(i, i) = stiffness_list[i];
	mKd(i, i) = damping_list[i];
    }

    for (const BodyNode *bn: bvh.skeleton->getBodyNodes())
	std::cout << bn->getName() << std::endl;

    std::vector<Eigen::VectorXd> frame;
    std::vector<std::vector<Eigen::Vector3d>> force;
    std::vector<Eigen::VectorXd> target;
    target.push_back(bvh.frame[startFrame]);
    for (size_t i = 0; i < 4; ++i)
	getTarget(bvh);
    std::vector<Eigen::VectorXd> tmp = getTarget(bvh);
    for (const Eigen::VectorXd &v: tmp)
    {
	for (size_t i = 0; i < groupNum; ++i)
	    target.push_back(v);
    }
    simulate(bvh, tmp, frame, force);
    std::vector<Eigen::VectorXd> ref;
    for (size_t i = 0; i < frame.size(); ++i)
	ref.push_back(bvh.frame[startFrame + i]);
    std::ofstream output;
    output.open(outputFileName);
    for (const Eigen::VectorXd &v: bvh.frameToEulerAngle(frame))
	output << v.transpose() << std::endl;
    output.close();
    output.open(outputFileName + ".raw.txt");
    for (const Eigen::VectorXd &v: frame)
	output << v.transpose() << std::endl;
    output.close();
    output.open(outputFileName + ".force.txt");
    output << force.size() << std::endl;
    for (const std::vector<Eigen::Vector3d> &list: force)
    {
	output << list.size() << std::endl;
	for (const Eigen::Vector3d &v: list)
	    output << v << std::endl;
    }
    output.close();
    output.open(outputFileName + ".target.txt");
    for (const Eigen::VectorXd &v: bvh.frameToEulerAngle(target))
	output << v.transpose() << std::endl;
    output.close();
    output.open(outputFileName + ".ref.txt");
    for (const Eigen::VectorXd &v: bvh.frameToEulerAngle(ref))
	output << v.transpose() << std::endl;
    output.close();
}
