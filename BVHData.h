#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cctype>
#include <Eigen/Core>
#include "dart/dart.hpp"

using namespace dart::dynamics;
using namespace dart::simulation;

class MyShapeNode
{
    public:
	std::shared_ptr<Shape> shape;
	Eigen::Isometry3d tf;

	MyShapeNode() {}
	MyShapeNode(std::shared_ptr<Shape> shape, Eigen::Isometry3d tf): shape(shape), tf(tf) {}
};

class BVHData
{
    public:
	SkeletonPtr skeleton;
	std::vector<Eigen::VectorXd> frame;
	std::vector<int> eulerAngleOrder;

	int loadBVH(const std::string& filename, const std::string& configFileName = "", const std::string& hingeJointFileName = "", double scale = 100.0);

	int parseBVH(std::istream& input);

	void setPositionAt(size_t i) const;

	size_t getChannelSize() const;

	std::vector<Eigen::VectorXd> frameToEulerAngle() const;

	std::vector<Eigen::VectorXd> frameToEulerAngle(std::vector<Eigen::VectorXd> targetFrame) const;

    private:
	size_t m_channelSize = 0;
	double scale;

	const std::vector<Eigen::Vector3d> axises = { Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY(), Eigen::Vector3d::UnitZ() };

	std::map<std::string, std::vector<MyShapeNode>> geometryConfig;

	std::map<std::string, int> hingeJoints;

	int loadGeometryConfig(std::istream &input);

	int parseRoot(std::istream& input);

	int parseIdentifier(std::istream& input, std::string& identifier);

	int parseToken(std::istream& input, const std::string& token);

	int parseOffset(std::istream& input, double& x, double& y, double& z);

	int parseChannels(std::istream& input, std::vector<std::string>& channels);

	int parseJoint(std::istream& input, const BodyNodePtr& parent, const std::string& name);

	int parseFrame(std::istream& input);

	BodyNodePtr addBody(const SkeletonPtr& skeleton, const BodyNodePtr& parent, const std::string& name, double x, double y, double z);

	void setGeometry(const BodyNodePtr& bn, double x, double y, double z, bool isEndSite = false);

	static void upper(std::string& s);

	static bool isIdentifier(const std::string& s);

	static std::vector<std::string> split(std::string s);
};
