#ifndef CONTROL_FRAGMENT_H
#define CONTROL_FRAGMENT_H

#include <vector>
#include <Eigen/Core>

class ControlFragment
{
    public:

	Eigen::MatrixXd m, sigma;
	Eigen::VectorXd a, tracked;
	Eigen::VectorXd vel; // pose = tracked

	std::vector<Eigen::Vector3d> endQ, endVel, endPos;
	Eigen::Vector3d endCOM, endCOMVel;
	Eigen::MatrixXd transformation;

	Eigen::VectorXd aftOffset;

};

#endif
