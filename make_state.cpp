#include <cmath>
#include <limits>
#include <cassert>
#include <Eigen/Geometry>
#include <dart/dart.hpp>
#include "utility.h"
#include "mass.h"

using namespace std;
using namespace Eigen;
using namespace dart::dynamics;

VectorXd makeState(const VectorXd &pose, const VectorXd &vel)
{
    SkeletonPtr skeleton = Utility::bvhs[omp_get_thread_num()].skeleton;
    skeleton->setPositions(pose);
    skeleton->setVelocities(vel);
    BodyNodePtr root = skeleton->getRootBodyNode();
    Isometry3d transform = root->getWorldTransform();
    Matrix3d rotation = transform.linear();
    Vector3d y = Vector3d::UnitY();
    Vector3d z = rotation.col(2);
    // we want to keep y axis verticall aligned and z axis aligned with the character's facing direction
    Vector3d x = y.cross(z);
    x.normalized();
    z = x.cross(y);
    z.normalized();
    transform.linear().col(0) = x;
    transform.linear().col(1) = y;
    transform.linear().col(2) = z;
    SimpleFrame frame(Frame::World());
    frame.setTransform(transform);

    // root transform in new frame
    transform = root->getTransform(&frame);
    rotation = transform.linear();
    AngleAxisd aa;
    aa.fromRotationMatrix(rotation);
    Vector3d axis = aa.axis();
    double angle = aa.angle();
    //assert(abs(axis.y()) <= numeric_limits<double>::epsilon());
    axis *= angle;
    Vector2d q0;
    q0 << axis.x(), axis.z();

    Vector3d translation = transform.translation();
    assert(abs(translation.x()) <= numeric_limits<double>::epsilon());
    assert(abs(translation.z()) <= numeric_limits<double>::epsilon());
    double h0 = translation.y();

    Vector6d c_and_c_dot = skeleton->getCOMSpatialVelocity(&frame);

    BodyNodePtr leftFoot = skeleton->getBodyNodes()[Utility::leftFootIndex];
    Vector3d dl = leftFoot->getCOM(&frame);

    BodyNodePtr rightFoot = skeleton->getBodyNodes()[Utility::rightFootIndex];
    Vector3d dr = rightFoot->getCOM(&frame);

    Vector3d l = Vector3d::Zero();
    for (BodyNodePtr bn: skeleton->getBodyNodes())
    {
	Vector3d am = bn->getAngularVelocity(&frame);
	am *= bn->getMass();
	l += am;
    }

    VectorXd state = VectorXd(18);
    state << q0, h0, c_and_c_dot, dl, dr, l;

    return state;
}
