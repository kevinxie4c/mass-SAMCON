#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <Eigen/Core>
#include <dart/dart.hpp>
#include "config_loader.h"
#include "BVHData.h"
#include "control_fragment.h"

namespace Utility
{
    extern BVHData mbvh;
    // do we really need to keep this?
    // make it OMP safe
    extern std::vector<BVHData> bvhs;
    extern Eigen::MatrixXd mKp, mKd;
    extern size_t ndof;

    void init();

    dart::dynamics::SkeletonPtr createFloor();

    // set position and velocity at index
    void setStateAt(size_t index, Eigen::VectorXd &pose, Eigen::VectorXd &vel);

    struct ErrorTerms
    {
	double err_p, err_r, err_e, err_b, err_zmp;
    };

    /* for improving SAMCON
       double costFunc(const dart::dynamics::SkeletonPtr, size_t index);
    */
    double costFunc(const dart::dynamics::SkeletonPtr, ControlFragment &cf, Eigen::Vector3d &zmp, ErrorTerms &et);

    template<typename T> std::vector<T> readListFrom(const std::string &filename)
    {
	std::ifstream input(filename);
	if (!input.good())
	    std::cerr << filename << ": not good!" << std::endl;
	std::vector<T> result;
	T d;
	while (input >> d)
	    result.push_back(d);
	return result;
	input.close();
    }

    bool fileGood(const std::string &filename);

    Eigen::VectorXd readVectorXdFrom(const std::string &filename);

    std::vector<std::string> split(std::string s);

    std::vector<Eigen::VectorXd> readVectorXdListFrom(const std::string &filename);

    Eigen::MatrixXd readMatrixXdFrom(const std::string &filename);

    bool dirExists(std::string name);
    void createDir(std::string name);

    extern std::vector<std::string> endEffectorName;
    extern size_t rootIndex;
    extern size_t waistIndex, leftHipIndex, rightHipIndex, leftKneeIndex, rightKneeIndex;
    extern size_t leftFootIndex, rightFootIndex;
    extern std::vector<size_t> endEffectorIndex;

    // for cost function
    extern std::vector<std::vector<Eigen::Vector3d>> target_q;
    extern std::vector<std::vector<Eigen::Vector3d>> target_omega;
    extern std::vector<std::vector<Eigen::Vector3d>> target_p;
    extern std::vector<Eigen::Vector3d> target_COM;
    extern std::vector<Eigen::Vector3d> target_COM_vel;
};

#endif
