#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <dart/dart.hpp>
#include "BVHData.h"
#include "utility.h"

using namespace std;
using namespace Eigen;
using namespace dart::dynamics;

SkeletonPtr skeleton;

VectorXd sum(const vector<VectorXd> &list)
{
    VectorXd v = VectorXd::Zero(list[0].size());
    for (auto u: list)
	v = skeleton->getPositionDifferences(v, -u);
    return v;
}

double var(const vector<VectorXd> &list)
{
    double res = 0;
    VectorXd m = sum(list) / list.size();
    for (auto u: list)
    {
	double n = skeleton->getPositionDifferences(u, m).norm();
	res += n * n;
    }
    res /= list.size();
    return res;
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
	cout << argv[0] << " bvh target reference" << endl;
	exit(0);
    }

    BVHData bvh;
    bvh.loadBVH(argv[1]);
    vector<VectorXd> tgt = bvh.frameToEulerAngle(Utility::readVectorXdListFrom(argv[2]));
    vector<VectorXd> ref = bvh.frameToEulerAngle(Utility::readVectorXdListFrom(argv[3]));
    skeleton = bvh.skeleton;
    vector<VectorXd> diff;
    for (size_t i = 0; i < tgt.size(); ++i)
	diff.push_back(skeleton->getPositionDifferences(tgt[i], ref[i]));
    cout << "NSR = " << var(diff) / var(ref) << endl;

    return 0;
}
