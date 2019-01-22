#include "BVHData.h"
#include <cmath>
#include <cfloat>
#include <utility>

int BVHData::loadBVH(const std::string& filename)
{
	std::ifstream input;
	frame.clear();
	eulerAngleOrder.clear();
	input.open(filename);
	int ret = parseBVH(input);
	input.close();
	skeleton->setSelfCollisionCheck(false);
	return ret;
}

int BVHData::parseBVH(std::istream& input)
{
	std::string tag;
	input >> tag;
	upper(tag);
	if (tag != "HIERARCHY")
	{
		std::cerr << "expected HIERARCHY but found " << tag << std::endl;
		return 0;
	}
	std::clog << "found HIERARCHY" << std::endl;
	parseRoot(input);
	parseFrame(input);
	return 1;
}


void BVHData::setPositionAt(size_t i) const
{
	skeleton->setPositions(frame[i]);
}

size_t BVHData::getChannelSize() const
{
	return m_channelSize;
}

std::vector<Eigen::VectorXd> BVHData::frameToEulerAngle() const
{
	return frameToEulerAngle(frame);
}

std::vector<Eigen::VectorXd> BVHData::frameToEulerAngle(std::vector<Eigen::VectorXd> targetFrame) const
{
	static double pi = acos(-1);
	std::vector<Eigen::VectorXd> list = targetFrame;
	std::vector<Eigen::VectorXd> result;
	for (Eigen::VectorXd &pose: list)
	{
		Eigen::VectorXd vec(getChannelSize());
		pose.segment(0, 3).swap(pose.segment(3, 3)); // swap position and rotation of the root
		pose.segment(0, 3) *= cm_per_m;
		vec.segment(0, 3) = pose.segment(0, 3);
		/*
		for (size_t i = 3; i < getChannelSize(); i += 3)
		{
			pose.segment(i, 3) = BallJoint::convertToRotation(pose.segment(i, 3)).eulerAngles(eulerAngleOrder[i], eulerAngleOrder[i + 1], eulerAngleOrder[i + 2]) / (2.0 * pi) * 360.0;
		}
		*/
		size_t index = 3;
		for (size_t i = 3; i < getChannelSize(); i += 3)
		{
			std::string jointName = skeleton->getJoints()[i / 3 - 1]->getName();
			bool isHinge = false;
			for (std::string &s: hingeJointList)
				if (jointName.find(s) != std::string::npos)
				{
					isHinge = true;
					break;
				}
			if (isHinge)
			{
				vec[i] = 0;
				vec[i + 1] = pose[index] / (2.0 * pi) * 360.0;
				vec[i + 2] = 0;
				index += 1;
			}
			else
			{
				vec.segment(i, 3) = BallJoint::convertToRotation(pose.segment(index, 3)).eulerAngles(eulerAngleOrder[i], eulerAngleOrder[i + 1], eulerAngleOrder[i + 2]) / (2.0 * pi) * 360.0;
				index += 3;
			}
		}
		result.push_back(vec);
	}
	return result;
}

int BVHData::parseRoot(std::istream& input)
{
	std::string tag, rootName;
	input >> tag;
	upper(tag);
	if (tag != "ROOT")
	{
		std::cerr << "expected ROOT but found " << tag << std::endl;
		return 0;
	}
	if (!parseIdentifier(input, rootName))
		return 0;
	std::clog << "found ROOT " << rootName << std::endl;
	skeleton = Skeleton::create("skeleton");
	if (!parseToken(input, "{"))
		return 0;
	parseJoint(input, nullptr, rootName);
	for (size_t i = 1; i < skeleton->getJoints().size(); ++i)
		skeleton->getJoints()[i]->setPositionLimitEnforced(true);
	static double pi = acos(-1);
	for (size_t i = 6; i < skeleton->getDofs().size(); ++i)
	{
		DegreeOfFreedom *dof = skeleton->getDofs()[i];
		dof->setPositionUpperLimit(pi);
		dof->setPositionLowerLimit(-pi);
	}
	return 1;
}

int BVHData::parseIdentifier(std::istream& input, std::string& identifier)
{
	input >> identifier;
	if (!isIdentifier(identifier))
	{
		std::cerr << "expected a identifier but found " << identifier << std::endl;
		return 0;
	}
	return 1;
}

int BVHData::parseToken(std::istream& input, const std::string& token)
{
	std::string tag;
	input >> tag;
	upper(tag);
	if (tag != token)
	{
		std::cerr << "expected a token \"" << token << "\" but found " << tag << std::endl;
		return 0;
	}
	return 1;
}

int BVHData::parseOffset(std::istream& input, double& x, double& y, double& z)
{
	if (input >> x >> y >> z)
		return 1;
	else
	{
		std::cerr << "invalid input" << std::endl;
		return 0;
	}
}

int BVHData::parseChannels(std::istream& input, std::vector<std::string>& channels)
{
	int n;
	if (input >> n)
	{
		channels.clear();
		for (int i = 0; i < n; ++i)
		{
			std::string channame;
			input >> channame;
			channels.push_back(channame);
		}
	}
	else
	{
		std::cerr << "invalid input" << std::endl;
		return 0;
	}
	return 1;
}
		

int BVHData::parseJoint(std::istream& input, const BodyNodePtr& parent, const std::string& name)
{
	std::string tag;
	bool offsetFound = false, channelsFound = false;
	BodyNodePtr bn = nullptr;
	while (input >> tag)
	{
		upper(tag);
		if (tag == "OFFSET")
		{
			if (offsetFound)
			{
				std::cerr << "redefined OFFSET" << std::endl;
				return 0;
			}
			offsetFound = true;
			std::clog << "found OFFSET" << std::endl;
			double x, y, z;
			if (!parseOffset(input, x, y, z)) 
			{
				std::cerr << "expected 3 doubles after OFFSET" << std::endl;
				return 0;
			}
			std::clog << "OFFSET " << x << " " << y << " " << z << std::endl;
			//bn = addBody(skeleton, parent, name, x, y, z);
			bn = addBody(skeleton, parent, name, x / cm_per_m, y / cm_per_m, z / cm_per_m);
		} 
		else if (tag == "CHANNELS")
		{
			if (channelsFound)
			{
				std::cerr << "redefined CHANNELS" << std::endl;
				return 0;
			}
			channelsFound = true;
			std::clog << "found CHANNELS" << std::endl;
			std::vector<std::string> channels;
			if (!parseChannels(input, channels))
			{
				std::cerr << "expected an int and a list of channel names after CHANNELS" << std::endl;
				return 0;
			}
			std::clog << "CHANNELS";
			for (std::string s: channels)
			{
				std::clog << " " << s;
				upper(s);
				if (s == "XROTATION")
					eulerAngleOrder.push_back(0);
				else if (s == "YROTATION")
					eulerAngleOrder.push_back(1);
				else if (s == "ZROTATION")
					eulerAngleOrder.push_back(2);
				else
					eulerAngleOrder.push_back(-1);
			}
			std::clog << std::endl;
			m_channelSize += channels.size();
		}
		else if (tag == "JOINT")
		{
			std::string jointName;
			if (!parseIdentifier(input, jointName))
				return 0;
			std::clog << "found JOINT " << jointName << std::endl;
			// TODO: make body
			if (!parseToken(input, "{"))
				return 0;
			if (bn == nullptr)
			{
				std::cerr << "bn == nullptr" << std::endl;
				return 1;
			}
			parseJoint(input, bn, jointName);
		}
		else if (tag == "END")
		{
			std::string t;
			input >> t;
			upper(t);
			if (t != "SITE")
			{
				std::cerr << "found END " << t << ", do you mean END SITE?" << std::endl;
				return 0;
			}
			std::clog << "found END SITE" << std::endl;
			if (!parseToken(input, "{"))
				return 0;
			if (!parseToken(input, "OFFSET"))
				return 0;
			double x, y, z;
			if (!parseOffset(input, x, y, z)) 
			{
				std::cerr << "expected 3 doubles after OFFSET" << std::endl;
				return 0;
			}
			std::clog << "OFFSET " << x << " " << y << " " << z << std::endl;
			if (bn == nullptr)
			{
				std::cerr << "bn == nullptr" << std::endl;
				return 0;
			}
			//setGeometry(bn, x, y, z, true);
			setGeometry(bn, x / cm_per_m, y / cm_per_m, z / cm_per_m, true);
			if (!parseToken(input, "}"))
				return 0;
		}
		else if (tag == "}")
		{
			if (!offsetFound)
			{
				std::cerr << "cannot found OFFSET before meeting \"}\"" << std::endl;
				return 0;
			}
			if (!channelsFound)
			{
				std::cerr << "cannot found CHANNELS before meeting \"}\"" << std::endl;
				return 0;
			}
			return 1;
		}
	}
	std::cerr << "cannot found \"}\"" << std::endl;
	return 0;
}

std::vector<std::string> BVHData::hingeJointList = { "Elbow", "Knee" };

int BVHData::parseFrame(std::istream& input)
{
	/*
	for (auto it: m_pos)
		std::cout << " " << it;
	std::cout << std::endl;
	*/
	static std::vector<Eigen::Vector3d> axises = { Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY(), Eigen::Vector3d::UnitZ() };
	if (!parseToken(input, "MOTION"))
		return 0;
	if (!parseToken(input, "FRAMES:"))
		return 0;
	size_t frameNum;
	if (input >> frameNum) // parse frame numbers
		std::cout << "frame number: " << frameNum << std::endl;
	else
	{
		std::cerr << "parse frame number failed" << std::endl;
		return 0;
	}
	if (!(parseToken(input, "FRAME") && parseToken(input, "TIME:")))
		return 0;
	double frameTime;
	if (input >> frameTime) // parse frame time
		std::cout << "frame time: " << frameTime << std::endl;
	else
	{
		std::cerr << "parse frame time failed" << std::endl;
		return 0;
	}
	for (size_t i = 0; i < frameNum; ++i)
	{
		Eigen::VectorXd t(getChannelSize());
		Eigen::VectorXd q(skeleton->getDofs().size());
		for (size_t j = 0; j < getChannelSize(); ++j)
		{
			double x;
			if (input >> x)
			{
				static double pi = acos(-1);
				if (j >= 3) // assume first 3 channels are position and the rest are rotation
					x = x / 360.0 * 2 * pi;	// degree to radius
				else
					x = x / cm_per_m; // for test
				t[j] = x;
			}
			else
			{
				std::cerr << "parse frame failed" << std::endl;
				return 0;
			}
		}
		size_t index = 3;
		for (size_t j = 3; j < getChannelSize(); j += 3)
		{
			bool isHinge = false;
			std::string jointName = skeleton->getJoints()[j / 3 - 1]->getName();
			for (std::string &s: hingeJointList)
				if (jointName.find(s) != std::string::npos)
				{
					isHinge = true;
					break;
				}
			if (isHinge)
			{
				q[index] = t[j + 1];
				index += 1;
			}
			else
			{
				Eigen::AngleAxisd r0 = Eigen::AngleAxisd(t[j], axises[eulerAngleOrder[j]]);
				Eigen::AngleAxisd r1 = Eigen::AngleAxisd(t[j + 1], axises[eulerAngleOrder[j + 1]]);
				Eigen::AngleAxisd r2 = Eigen::AngleAxisd(t[j + 2], axises[eulerAngleOrder[j + 2]]);
				Eigen::Matrix3d rotation = (r0 * r1 * r2).toRotationMatrix();
				Eigen::Vector3d positions = BallJoint::convertToPositions(rotation);
				q[index] = positions[0];
				q[index + 1] = positions[1];
				q[index + 2] = positions[2];
				index += 3;
			}
		}
		for (size_t j = 0; j < 3; ++j)
		    q[j] = t[j];
		// root node: rotation, translation
		for (size_t j = 0; j < 3; ++j)
			std::swap(q[j], q[j + 3]);
		frame.push_back(q);
	}
	return 1;
}

BodyNodePtr BVHData::addBody(const SkeletonPtr& skeleton, const BodyNodePtr& parent, const std::string& name, double x, double y, double z)
{
	BallJoint::Properties ballJointProperties;
	ballJointProperties.mName = name + "_joint";
	ballJointProperties.mT_ParentBodyToJoint.translation() = Eigen::Vector3d(x, y, z);
	
	FreeJoint::Properties freeJointProperties;
	freeJointProperties.mName = name + "_joint";
	freeJointProperties.mT_ParentBodyToJoint.translation() = Eigen::Vector3d(x, y, z);

	RevoluteJoint::Properties revoluteJointProperties;
	revoluteJointProperties.mName = name + "_joint";
	revoluteJointProperties.mT_ParentBodyToJoint.translation() = Eigen::Vector3d(x, y, z);
	revoluteJointProperties.mAxis = Eigen::Vector3d::UnitX();

	BodyNodePtr bn;
	if (parent == nullptr) 
		bn = skeleton->createJointAndBodyNodePair<FreeJoint>(parent, freeJointProperties, BodyNode::AspectProperties(name)).second;	// Assume that root is a free joint
	else
	{
		bool isHinge = false;
		for (std::string &s: BVHData::hingeJointList)
			if (name.find(s) != std::string::npos)
			{
				isHinge = true;
				break;
			}
		if (isHinge)
			bn = skeleton->createJointAndBodyNodePair<RevoluteJoint>(parent, revoluteJointProperties, BodyNode::AspectProperties(name)).second;
		else
			bn = skeleton->createJointAndBodyNodePair<BallJoint>(parent, ballJointProperties, BodyNode::AspectProperties(name)).second;
	}
	bn->setFrictionCoeff(1e20);
	bn->setRestitutionCoeff(0);

	// Make a shape for the Joint
	const double& R = 0.02;
	std::shared_ptr<EllipsoidShape> ball(
		new EllipsoidShape(sqrt(2) * Eigen::Vector3d(R, R, R)));
	auto shapeNode = bn->createShapeNodeWith<VisualAspect>(ball);
	shapeNode->getVisualAspect()->setColor(dart::Color::Blue());

	// TODO: set geometry
	if (parent != nullptr)
		setGeometry(parent, x, y, z);
	return bn;
}

void BVHData::setGeometry(const BodyNodePtr& bn, double x, double y, double z, bool isEndSite)
{
	double w = 0.02;
	double scale = 1.0;
	if (isEndSite) scale = 1.5;
	std::shared_ptr<BoxShape> box(new BoxShape(Eigen::Vector3d(scale * sqrt(x * x + y * y + z * z), w, w)));
	auto shapeNode = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(box);
	// TODO
	shapeNode->getVisualAspect()->setColor(dart::Color::Blue());
	Eigen::Isometry3d tf;
	Eigen::Vector3d t(x / 2.0, y / 2.0, z / 2.0);
	Eigen::Vector3d xr(x, y, z);
	xr.normalize();
	Eigen::Vector3d yr;
	if (std::abs(x) < DBL_EPSILON && std::abs(z) < DBL_EPSILON)
		yr << 0, 0, 1;
	else
		yr << 0, 1, 0;
	Eigen::Vector3d zr = xr.cross(yr);
	zr.normalize();
	yr = zr.cross(xr);
	yr.normalize();
	tf.matrix() << xr, yr, zr, t, 0, 0, 0, 1;
	shapeNode->setRelativeTransform(tf);
	bn->setLocalCOM(t);
}

void BVHData::upper(std::string& s)
{
	for (auto& c: s) c = toupper(c);
}

bool BVHData::isIdentifier(const std::string& s)
{
	for (auto& c: s) if (!(('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z') || c == '_')) return false;
	return true;
}
