#include <iostream>
#include <algorithm>
#include <cmath>
#include <dart/collision/ode/OdeCollisionDetector.hpp>
#include <sys/stat.h>
#include "utility.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#endif

// make setStateAt OMP safe
static std::vector<dart::dynamics::SkeletonPtr> skeletons, t_skeletons; // for velocity computation
/*
static dart::dynamics::SkeletonPtr m_floor;
static dart::collision::CollisionOption option;
static dart::collision::CollisionResult result;
static dart::collision::CollisionDetectorPtr collisionEngine = dart::collision::OdeCollisionDetector::create();
*/

namespace Utility
{
    BVHData mbvh;
    std::vector<BVHData> bvhs;

    Eigen::MatrixXd mKp, mKd;
    size_t ndof;

    // TODO: change name: xxxName -> xxxNameList
    std::vector<std::string> endEffectorName = { "Wrist", "Head", "Ankle", "Foot" };
    std::vector<std::string> leftFootNameList = { "LeftAngle", "LeftFoot" };
    std::vector<std::string> rightFootNameList = { "RightAnkle", "RightFoot" };
    std::vector<std::string> waistName = { "Abdomen" };
    std::vector<std::string> leftHipNameList = { "LeftHip" };
    std::vector<std::string> rightHipNameList = { "RightHip" };
    std::vector<std::string> leftKneeName = { "LeftKnee" };
    std::vector<std::string> rightKneeName = { "RightKnee" };
    std::string leftFootName;
    std::string rightFootName;
    size_t rootIndex = 0;
    size_t leftFootIndex, rightFootIndex;
    size_t leftHipDofIdx, rightHipDofIdx;
    size_t waistIndex, leftHipIndex, rightHipIndex, leftKneeIndex, rightKneeIndex;
    std::vector<size_t> endEffectorIndex;

    // for cost function
    std::vector<std::vector<Eigen::Vector3d>> target_q;
    std::vector<std::vector<Eigen::Vector3d>> target_omega;
    std::vector<std::vector<Eigen::Vector3d>> target_p;
    std::vector<Eigen::Vector3d> target_COM;
    std::vector<Eigen::Vector3d> target_COM_vel;
};

void Utility::init()
{
    dart::constraint::ContactConstraint::setErrorReductionParameter(Config::ERP);
    dart::constraint::ContactConstraint::setConstraintForceMixing(Config::CFM);

    // inititialize name lists
    if (Config::endEffectorFileName != "")
	endEffectorName= readListFrom<std::string>(Config::endEffectorFileName);
    if (Config::leftFootFileName != "")
	leftFootNameList = readListFrom<std::string>(Config::leftFootFileName);
    if (Config::rightFootFileName != "")
	rightFootNameList = readListFrom<std::string>(Config::rightFootFileName);
    if (Config::waistFileName != "")
	waistName = readListFrom<std::string>(Config::waistFileName);
    if (Config::leftHipFileName != "")
	leftHipNameList = readListFrom<std::string>(Config::leftHipFileName);
    if (Config::rightHipFileName != "")
	rightHipNameList = readListFrom<std::string>(Config::rightHipFileName);
    if (Config::leftKneeFileName != "")
	leftKneeName = readListFrom<std::string>(Config::leftKneeFileName);
    if (Config::rightKneeFileName != "")
	rightKneeName = readListFrom<std::string>(Config::rightKneeFileName);


    mbvh.loadBVH(Config::bvhFileName, Config::geometryConfigFileName, Config::hingeJointListFileName, Config::scale);
    ndof = mbvh.skeleton->getDofs().size();
    std::cout << "ndof: " << ndof << std::endl;
    std::vector<double> massList = readListFrom<double>(Config::massFileName);
    std::vector<dart::dynamics::BodyNode*> bns = mbvh.skeleton->getBodyNodes();
    std::cout << "massList size: " << massList.size() << std::endl;
    std::cout << "bns size: " << bns.size() << std::endl;
    if (massList.size() != bns.size())
	std::cerr << "massList.size() != bns.size()" << std::endl;
    
    // set mass
    for (size_t i = 0; i < massList.size(); ++i)
    {
	std::cout << bns[i]->getName() << " " << massList[i] << std::endl;
	/*
	std::cout << bns[i]->getShapeNodes().size() << std::endl;
	ShapeNode *shapeNode = bns[i]->getShapeNode(1);
	ShapePtr shape = shapeNode->getShape();
	std::cout << shape->getType() << std::endl;
	std::cout << shape->getVolume() << std::endl;
	Inertia inertia;
	inertia.setLocalCOM(bns[i]->getLocalCOM());
	inertia.setMass(massList[i]);
	inertia.setMoment(shape->computeInertia(massList[i]));
	std::cout << inertia.getLocalCOM() << std::endl;
	std::cout << inertia.getMoment() << std::endl;
	bns[i]->setInertia(inertia);
	*/
	std::cout << bns[i]->getSpatialInertia() << std::endl;
	bns[i]->setMass(massList[i]);
    }
    std::vector<double> stiffness = readListFrom<double>(Config::stiffnessFileName);
    std::vector<double> damping = readListFrom<double>(Config::dampingFileName);
    
    mKp = Eigen::MatrixXd::Zero(ndof, ndof);
    mKd = Eigen::MatrixXd::Zero(ndof, ndof);
    for (size_t i = 0; i < ndof; ++i)
    {
	mKp(i, i) = stiffness[i];
	mKd(i, i) = damping[i];
    }

    for (size_t i = 0; i < omp_get_max_threads(); ++i)
    {
	skeletons.push_back(mbvh.skeleton->clone());
	t_skeletons.push_back(mbvh.skeleton->clone());
    }

    // init cost
    // not using this one now
    for (size_t index = 0; index < mbvh.frame.size(); ++index)
    {
	Eigen::VectorXd pose;
	Eigen::VectorXd vel;
	if (index < mbvh.frame.size() - 1)
	    setStateAt(index, pose, vel);
	else
	{
	    // use velocity at index - 1 for last frame
	    setStateAt(index - 1, pose, vel);
	    pose = mbvh.frame.back();
	}
	mbvh.skeleton->setPositions(pose);
	mbvh.skeleton->setVelocities(vel);
	std::vector<Eigen::Vector3d> q;
	std::vector<Eigen::Vector3d> omega;
	for (size_t i = 0; i < mbvh.skeleton->getJoints().size(); ++i)
	{
	    // root has 6 DOFs, so?
	    Joint *joint = mbvh.skeleton->getJoints()[i];
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
	target_q.push_back(q);
	target_omega.push_back(omega);
	std::vector<Eigen::Vector3d> p;
	for (BodyNode const *bodyNode: mbvh.skeleton->getBodyNodes())
	    p.push_back(bodyNode->getCOM());
	target_p.push_back(p);
	target_COM.push_back(mbvh.skeleton->getCOM());
	target_COM_vel.push_back(mbvh.skeleton->getCOMLinearVelocity());
    }

    // endEffectorIndex
    for (size_t i = 0; i < mbvh.skeleton->getBodyNodes().size(); ++i)
    {
	std::cout << mbvh.skeleton->getBodyNodes()[i]->getName() << std::endl;
	for (const std::string& str: endEffectorName)
	    if (mbvh.skeleton->getBodyNodes()[i]->getName().find(str) != std::string::npos)
	    {
		endEffectorIndex.push_back(i);
		break;
	    }
    }
    std::cout << "endEffectorIndex:";
    for (size_t &i: endEffectorIndex)
	std::cout << " " << i;
    std::cout << std::endl;

    bool found = false;
    for (size_t i = 0; i < mbvh.skeleton->getBodyNodes().size() && !found; ++i)
    {
	for (const std::string& str: leftFootNameList)
	    if (mbvh.skeleton->getBodyNodes()[i]->getName().find(str) != std::string::npos)
	    {
		leftFootIndex = i;
		leftFootName = mbvh.skeleton->getBodyNodes()[i]->getName();
		found = true;
		break;
	    }
    }
    found = false;
    for (size_t i = 0; i < mbvh.skeleton->getBodyNodes().size() && !found; ++i)
    {
	for (const std::string& str: rightFootNameList)
	    if (mbvh.skeleton->getBodyNodes()[i]->getName().find(str) != std::string::npos)
	    {
		rightFootIndex = i;
		rightFootName = mbvh.skeleton->getBodyNodes()[i]->getName();
		found = true;
		break;
	    }
    }
    found = false;
    for (size_t i = 0; i < mbvh.skeleton->getDofs().size() && !found; ++i)
    {
	for (const std::string& str: leftHipNameList)
	    if (mbvh.skeleton->getDofs()[i]->getName().find(str) != std::string::npos)
	    {
		leftHipDofIdx = i;
		found = true;
	    }
    }
    found = false;
    for (size_t i = 0; i < mbvh.skeleton->getDofs().size() && !found; ++i)
    {
	for (const std::string& str: rightHipNameList)
	    if (mbvh.skeleton->getDofs()[i]->getName().find(str) != std::string::npos)
	    {
		rightHipDofIdx = i;
		found = true;
	    }
    }
    std::cout << "leftFootIndex = " << leftFootIndex << std::endl;
    std::cout << "rightFootIndex = " << rightFootIndex << std::endl;
    std::cout << "leftFootDofIdx = " << leftHipDofIdx << std::endl;
    std::cout << "rightFootDofIdx = " << rightHipDofIdx << std::endl;

    found = false;
    for (size_t i = 0; i < mbvh.skeleton->getJoints().size() && !found; ++i)
    {
	for (const std::string& str: waistName)
	    if (mbvh.skeleton->getJoint(i)->getName().find(str) != std::string::npos)
	    {
		waistIndex = i;
		found = true;
		break;
	    }
    }
    found = false;
    for (size_t i = 0; i < mbvh.skeleton->getJoints().size() && !found; ++i)
    {
	for (const std::string& str: leftHipNameList)
	    if (mbvh.skeleton->getJoint(i)->getName().find(str) != std::string::npos)
	    {
		leftHipIndex = i;
		found = true;
		break;
	    }
    }
    found = false;
    for (size_t i = 0; i < mbvh.skeleton->getJoints().size() && !found; ++i)
    {
	for (const std::string& str: rightHipNameList)
	    if (mbvh.skeleton->getJoint(i)->getName().find(str) != std::string::npos)
	    {
		rightHipIndex = i;
		found = true;
		break;
	    }
    }
    found = false;
    for (size_t i = 0; i < mbvh.skeleton->getJoints().size() && !found; ++i)
    {
	for (const std::string& str: leftKneeName)
	    if (mbvh.skeleton->getJoint(i)->getName().find(str) != std::string::npos)
	    {
		leftKneeIndex = i;
		found = true;
		break;
	    }
    }
    found = false;
    for (size_t i = 0; i < mbvh.skeleton->getJoints().size() && !found; ++i)
    {
	for (const std::string& str: rightKneeName)
	    if (mbvh.skeleton->getJoint(i)->getName().find(str) != std::string::npos)
	    {
		rightKneeIndex = i;
		found = true;
		break;
	    }
    }
    std::cout << "waistIndex = " << waistIndex << std::endl;
    std::cout << "leftHipIndex = " << leftHipIndex << std::endl;
    std::cout << "rightHipIndex = " << rightHipIndex << std::endl;
    std::cout << "leftKneeIndex = " << leftKneeIndex << std::endl;
    std::cout << "rightKneeIndex = " << rightKneeIndex << std::endl;

    for (size_t i = 0; i < omp_get_max_threads(); ++i)
	bvhs.push_back(mbvh);
    //m_floor = createFloor();
}

dart::dynamics::SkeletonPtr Utility::createFloor()
{
    using namespace dart::dynamics;
    SkeletonPtr floor = Skeleton::create("floor");

    // Give the floor a body
    BodyNodePtr body = floor->createJointAndBodyNodePair<WeldJoint>(nullptr).second;
    body->setFrictionCoeff(Config::frictionCoeff);
    body->setRestitutionCoeff(Config::restitutionCoeff);
    body->setName("floor");

    // Give the body a shape
    double floor_width = 1e8;
    double floor_height = 0.01;
    std::shared_ptr<BoxShape> box(
	    new BoxShape(Eigen::Vector3d(floor_width, floor_height, floor_width)));
    auto shapeNode = body->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(box);
    shapeNode->getVisualAspect()->setColor(dart::Color::Black());

    // Put the body into position
    Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
    tf.translation() = Eigen::Vector3d(0.0, Config::groundOffset - floor_height / 2, 0.0);
    body->getParentJoint()->setTransformFromParentBodyNode(tf);
    return floor;
}

// set position and velocity at index
// OMP safe (maybe)
void Utility::setStateAt(size_t index, Eigen::VectorXd &pose, Eigen::VectorXd &vel)
{
    SkeletonPtr skeleton = skeletons[omp_get_thread_num()], t_skeleton = t_skeletons[omp_get_thread_num()];
    skeleton->setPositions(mbvh.frame[index]);
    t_skeleton->setPositions(mbvh.frame[index + 1]);
    // also need to set velocities
    // what's the rule when we use rotation vector?
    // assume that the first one is the root joint
    {
	// root joint
	/*
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
	*/
	
	Joint *joint = skeleton->getJoint(0);
	Joint *t_joint = t_skeleton->getJoint(0);
	Eigen::VectorXd q1 = joint->getPositions();
	Eigen::VectorXd q2 = t_joint->getPositions();
	Eigen::VectorXd v = joint->getPositionDifferences(q2, q1) / Config::frameTime;
	if (Config::zeroInitialAngularVelocities)
	    v.head(3) = Eigen::Vector3d::Zero();
	joint->setVelocities(v);
    }
    // other joints
    for (size_t i = 1; i < skeleton->getJoints().size(); ++i)
    {
	/*
	Joint *joint = skeleton->getJoint(i);
	Eigen::Vector3d omega;
	if (joint->getPositions().rows() == 3)
	    omega = joint->getPositions();
	else
	    omega = Eigen::Vector3d(joint->getPosition(0), 0, 0);
	Joint *t_joint = t_skeleton->getJoint(i);
	Eigen::Vector3d t_omega;
	if (t_joint->getPositions().rows() == 3)
	    t_omega = t_joint->getPositions();
	else
	    t_omega = Eigen::Vector3d(t_joint->getPosition(0), 0, 0);
	Eigen::AngleAxisd r(omega.norm(), omega.normalized());
	Eigen::AngleAxisd t_r(t_omega.norm(), t_omega.normalized());
	Eigen::AngleAxisd diff_r(t_r * r.inverse());
	if (joint->getPositions().rows() == 3)
	    joint->setVelocities(diff_r.axis() * (diff_r.angle() / frameTime));
	else
	    joint->setVelocity(0, ((t_joint->getPosition(0) - joint->getPosition(0)) / frameTime));
	*/
	Joint *joint = skeleton->getJoint(i);
	Joint *t_joint = t_skeleton->getJoint(i);
	Eigen::VectorXd q1 = joint->getPositions();
	Eigen::VectorXd q2 = t_joint->getPositions();
	Eigen::VectorXd v = joint->getPositionDifferences(q2, q1) / Config::frameTime;
	joint->setVelocities(v);
    }
    pose = skeleton->getPositions();
    vel = skeleton->getVelocities();
}

/* for improving SAMCON
double Utility::costFunc(const SkeletonPtr skeleton, size_t index)
{
    std::vector<dart::dynamics::Joint*> joints = skeleton->getJoints();
    static size_t n = joints.size();
    double err_p = 0;
    std::vector<Eigen::Vector3d> qlist = target_q[index];
    std::vector<Eigen::Vector3d> omegalist = target_omega[index];
    for (size_t i = 0; i < joints.size(); ++i)
    {
	const Joint *joint = joints[i];
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
	    if (joint->getPositions().rows() == 3)
	    {
		q = joint->getPositions();
		omega= joint->getVelocities();
	    }
	    else
	    {
		q = Eigen::Vector3d(joint->getPosition(0), 0, 0);
		omega = Eigen::Vector3d(joint->getVelocity(0), 0, 0);
	    }
	}
	Eigen::Quaterniond quat1(Eigen::AngleAxisd(q.norm(), q.normalized()));
	Eigen::Quaterniond quat2(Eigen::AngleAxisd(qr.norm(), qr.normalized()));
	err_p += fabs(quat1.angularDistance(quat2)) + (omega - omega_r).norm() * 0.1;
    }
    err_p /= n;

    double err_r = 0;
    {
	const Joint *joint = joints[rootIndex];
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

    static size_t k = endEffectorIndex.size();
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

    return Config::wp * err_p + Config::wr * err_r + Config::we * err_e + Config::wb * err_b;
}
*/

double Utility::costFunc(const SkeletonPtr skeleton, ControlFragment &cf, Eigen::Vector3d &zmp, ErrorTerms &et)
{
    std::vector<dart::dynamics::Joint*> joints = skeleton->getJoints();
    static size_t n = joints.size();
    double err_p = 0;
    std::vector<Eigen::Vector3d> qlist = cf.endQ;
    std::vector<Eigen::Vector3d> omegalist = cf.endVel;
    for (size_t i = 0; i < joints.size(); ++i)
    {
	const Joint *joint = joints[i];
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
	    if (joint->getPositions().rows() == 3)
	    {
		q = joint->getPositions();
		omega= joint->getVelocities();
	    }
	    else
	    {
		q = Eigen::Vector3d(joint->getPosition(0), 0, 0);
		omega = Eigen::Vector3d(joint->getVelocity(0), 0, 0);
	    }
	}
	Eigen::Quaterniond quat1(Eigen::AngleAxisd(q.norm(), q.normalized()));
	Eigen::Quaterniond quat2(Eigen::AngleAxisd(qr.norm(), qr.normalized()));
	err_p += fabs(quat1.angularDistance(quat2)) + (omega - omega_r).norm() * 0.1;
    }
    err_p /= n;

    double err_r = 0;
    {
	const Joint *joint = joints[rootIndex];
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

    static size_t k = endEffectorIndex.size();
    double err_e = 0;
    std::vector<Eigen::Vector3d> plist = cf.endPos;
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
    Eigen::Vector3d COMr = cf.endCOM;
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
    err_b += (skeleton->getCOMLinearVelocity() - cf.endCOMVel).norm() * 0.1;

    /*
    auto leftFootGroup = collisionEngine->createCollisionGroup(skeleton->getBodyNode(leftFootIndex));
    auto rightFootGroup = collisionEngine->createCollisionGroup(skeleton->getBodyNode(rightFootIndex));
    auto floorGroup = collisionEngine->createCollisionGroup(m_floor.get());
    bool leftCollision = floorGroup->collide(leftFootGroup.get(), option, &result);
    bool rightCollision = floorGroup->collide(rightFootGroup.get(), option, &result);
    */

    // TODO: confirm the following computation of dL is correct
    Eigen::Vector3d dL = Eigen::Vector3d::Zero(); // rate of angular momentum at center of mass
    Eigen::Vector3d com = skeleton->getCOM();
    for (const BodyNode *bn: skeleton->getBodyNodes())
    {
	dL += (bn->getCOM() - com).cross(bn->getMass() * bn->getCOMLinearAcceleration());

	Eigen::Matrix3d R = bn->getTransform().rotation();
	Eigen::Matrix3d I = bn->getInertia().getMoment(); // moment of inertia in local frame
	Eigen::Vector3d omiga = bn->getCOMSpatialVelocity().head(3); // angular velocity in local frame
	Eigen::Vector3d dOmiga = bn->getCOMSpatialAcceleration().head(3); // angular acceleration in local frame

	dL += R * (I * dOmiga + omiga.cross(I * omiga));
    }

    double m = skeleton->getMass();
    Eigen::Vector3d g(0, Config::gravity, 0);
    Eigen::Vector3d F_gi = m * (g - skeleton->getCOMLinearAcceleration());
    Eigen::Vector3d P(com.x(), Config::groundOffset, com.z());
    Eigen::Vector3d vPG = com - P;
    Eigen::Vector3d M_gi_P = vPG.cross(m * g) - vPG.cross(m * skeleton->getCOMLinearAcceleration()) - dL;
    Eigen::Vector3d vn(0, 1, 0);
    Eigen::Vector3d vPZ = vn.cross(M_gi_P) / F_gi.dot(vn);
    //Eigen::Vector3d zmp = P + vPZ;
    zmp = P + vPZ;
    Eigen::Vector3d pL = skeleton->getBodyNode(leftFootIndex)->getCOM();
    Eigen::Vector3d pR = skeleton->getBodyNode(rightFootIndex)->getCOM();
    pL.y() = pR.y() = Config::groundOffset;
    double d = std::min(
	    std::min((zmp - pL).norm(), (zmp - pR).norm()),
	    std::abs((pR - pL).normalized().dot(zmp - pL))
	    );

    double err_zmp = d > 0.2 ? d - 0.2  : 0;
    et.err_p = err_p;
    et.err_r = err_r;
    et.err_e = err_e;
    et.err_b = err_b;
    et.err_zmp = err_zmp;
    return Config::wp * err_p + Config::wr * err_r + Config::we * err_e + Config::wb * err_b + Config::w_zmp * err_zmp;
}

bool Utility::fileGood(const std::string &filename)
{
    std::ifstream file(filename);
    return file.good();
}

Eigen::VectorXd Utility::readVectorXdFrom(const std::string &filename)
{
    std::vector<double> list = Utility::readListFrom<double>(filename);
    Eigen::VectorXd result(list.size());
    for (size_t i = 0; i < list.size(); ++i)
	result[i] = list[i];
    return result;
}

std::vector<std::string> Utility::split(std::string s)
{
    std::string delimiter = " ";
    std::vector<std::string> result;
    size_t pos;
    while ((pos = s.find(delimiter)) != std::string::npos)
    {
	std::string t = s.substr(0, pos);
	if (!t.empty())
	    result.push_back(t);
	s.erase(0, pos + delimiter.length());
    }
    result.push_back(s);
    return result;
}

std::vector<Eigen::VectorXd> Utility::readVectorXdListFrom(const std::string &filename)
{
    std::vector<Eigen::VectorXd> result;
    std::ifstream input(filename);
    if (!input.good())
	std::cerr << filename << ": not good!" << std::endl;
    std::string line;
    while (std::getline(input, line))
    {
	std::vector<std::string> list = split(line);
	Eigen::VectorXd vec(list.size());
	for (size_t i = 0; i < list.size(); ++i)
	    if (std::isspace(list[i][0]))
		break;
	    else 
		vec[i] = std::stod(list[i]);
	result.push_back(vec);
    }
    input.close();
    return result;
}

Eigen::MatrixXd Utility::readMatrixXdFrom(const std::string &filename)
{
    std::vector<std::vector<double>> list2d;
    std::ifstream input(filename);
    if (!input.good())
	std::cerr << filename << ": not good!" << std::endl;
    std::string line;
    size_t row = 0, col = 0;
    while (std::getline(input, line))
    {
	++row;
	std::vector<std::string> list = split(line);
	if (col == 0)
	    col = list.size();
	else
	{
	    if (col != list.size())
		std::cerr << "doesn't look like a matrix" << std::endl;
	}
	list2d.push_back(std::vector<double>());
	for (size_t i = 0; i < list.size(); ++i)
	    if (std::isspace(list[i][0]))
		break;
	    else 
		list2d.back().push_back(std::stod(list[i]));
    }
    input.close();
    Eigen::MatrixXd result(row, col);
	for (size_t i = 0; i < row; ++i)
	    for (size_t j = 0; j < col; ++j)
		result(i, j) = list2d[i][j];

    return result;
}

bool Utility::dirExists(std::string name)
{
    struct stat buf;
    if (stat(name.c_str(), &buf) == 0 && S_ISDIR(buf.st_mode))
	return true;
    return false;
}

void Utility::createDir(std::string name)
{
    mkdir(name.c_str(), 0755);
}
    
