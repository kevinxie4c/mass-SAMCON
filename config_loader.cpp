#include "config_loader.h"
#include <iostream>
#include <fstream>

namespace Config
{
	// Simulation
	double frictionCoeff = 0.5;
	double restitutionCoeff = 0.0;
	double gravity = -9.8;
	double ERP = 0.2;
	double CFM = 0.01;
	double timeStep = 0.0001;
	double frameTime = 0.008333;
	size_t stepPerFrame;

	// cost function
	double wp= 5, wr = 3, we = 30, wb = 10;

	// file name
	std::string bvhFileName;
	std::string geometryConfigFileName;
	std::string hingeJointListFileName;
	std::string freeNodesListFileName;
	std::string outputFileName = "result.txt";
	std::string basedir;
	std::string scalesFileName;
	std::string covFileName;
	std::string stiffnessFileName;
	std::string dampingFileName;
	std::string massFileName;
	std::string initMeanFileName;
	std::string endEffectorFileName;
	std::string leftFootFileName;
	std::string rightFootFileName;
	std::string waistFileName;
	std::string leftHipFileName;
	std::string rightHipFileName;
	std::string leftKneeFileName;
	std::string rightKneeFileName;

	// algorithm
	double scale = 100.0;
	size_t startFrame = 0;
	size_t endFrame = 1e8;
	double groundOffset = -0.03;
	size_t sampleNum = 1200, saveNum = 400;
	double percentageDiscard = 0.4;
	int convexExp = 6;
	size_t groupNum = 5;
	size_t slidingWindow = 50;
	size_t updateWindow = 50;
	size_t trialMin = 5, trialMax = 20;
	size_t notImproveMax = 5;
	double goodEnough = 5.0, failThreshold = 20.0;
	bool stablePD = true;
	bool zeroInitialVelocities = true;
	bool zeroInitialAngularVelocities = true;
	bool useSampleNumAsLambda = false;
	bool useAFT = false;
	double initSigma = 0.1;
	size_t maxFailure = 1;
	size_t rank = 10;
	size_t dRank = 5;
	double scaleMassMatrix;
	bool autoScaleMassMatrix = false;
	bool useEigenvalueScale = false;
	bool useInverseForce = false;
	bool useEigenvector = false;
	bool useEigenvalue = false;

	// others
	bool showWindow = false;
	bool onlyLogAndFinal = false;
	size_t loopNum = 50;
	bool generateSamplesFile = false;

};

void Config::load(const std::string &filename)
{
    std::ifstream file(filename);
    std::string parameter;
    std::string value;
    while (file >> parameter >> value)
    {
	std::cout << parameter << " = " << value << std::endl;
	setParameter(parameter, value);
    }
    file.close();
    stepPerFrame = frameTime / timeStep;
}

void Config::setParameter(const std::string &parameter, const std::string &value)
{
    if (parameter == "bvhFileName")
	bvhFileName = value;
    else if (parameter == "geometryConfigFileName")
	geometryConfigFileName = value;
    else if (parameter == "hingeJointListFileName")
	hingeJointListFileName = value;
    else if (parameter == "freeNodesListFileName")
	freeNodesListFileName = value;
    else if (parameter == "outputFileName")
	outputFileName = value;
    else if (parameter == "scalesFileName")
	scalesFileName = value;
    else if (parameter == "initMeanFileName")
	initMeanFileName = value;
    else if (parameter == "basedir")
	basedir = value;
    else if (parameter == "scale")
	scale = std::stod(value);
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
    else if (parameter == "percentageDiscard")
	percentageDiscard = std::stod(value);
    else if (parameter == "convexExp")
	convexExp = std::stoi(value);
    else if (parameter == "groupNum")
	groupNum = std::stoi(value);
    else if (parameter == "frameTime")
	frameTime = std::stod(value);
    else if (parameter == "timeStep")
	timeStep = std::stod(value);
    else if (parameter == "maxFailure")
	maxFailure = std::stoi(value);
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
    else if (parameter == "covFileName")
	covFileName = value;
    else if (parameter == "stiffnessFileName")
	stiffnessFileName = value;
    else if (parameter == "dampingFileName")
	dampingFileName = value;
    else if (parameter == "massFileName")
	massFileName = value;
    else if (parameter == "gravity")
	gravity = std::stod(value);
    else if (parameter == "stablePD")
	stablePD = std::stoi(value);
    else if (parameter == "zeroInitialVelocities")
	zeroInitialVelocities = std::stoi(value);
    else if (parameter == "zeroInitialAngularVelocities")
	zeroInitialAngularVelocities = std::stoi(value);
    else if (parameter == "showWindow")
	showWindow = std::stoi(value);
    else if (parameter == "onlyLogAndFinal")
	onlyLogAndFinal = std::stoi(value);
    else if (parameter == "ERP")
	ERP = std::stod(value);
    else if (parameter == "CFM")
	CFM = std::stod(value);
    else if (parameter == "slidingWindow")
	slidingWindow = std::stoi(value);
    else if (parameter == "updateWindow")
	updateWindow = std::stoi(value);
    else if (parameter == "useSampleNumAsLambda")
	useSampleNumAsLambda = std::stoi(value);
    else if (parameter == "useAFT")
	useAFT = std::stoi(value);
    else if (parameter == "initSigma")
	initSigma = std::stod(value);
    else if (parameter == "failThreshold")
	failThreshold = std::stod(value);
    else if (parameter == "goodEnough")
	goodEnough = std::stod(value);
    else if (parameter == "rank")
	rank = std::stoi(value);
    else if (parameter == "dRank")
	dRank = std::stoi(value);
    else if (parameter == "scaleMassMatrix")
	scaleMassMatrix = std::stod(value);
    else if (parameter == "autoScaleMassMatrix")
	autoScaleMassMatrix = std::stod(value);
    else if (parameter == "useEigenvalueScale")
	useEigenvalueScale = std::stoi(value);
    else if (parameter == "useInverseForce")
	useInverseForce = std::stoi(value);
    else if (parameter == "useEigenvector")
	useEigenvector = std::stoi(value);
    else if (parameter == "useEigenvalue")
	useEigenvalue = std::stoi(value);
    else if (parameter == "trialMax")
	trialMax = std::stoi(value);
    else if (parameter == "trialMin")
	trialMin = std::stoi(value);
    else if (parameter == "notImproveMax")
	notImproveMax = std::stoi(value);
    else if (parameter == "loopNum")
	loopNum = std::stoi(value);
    else if (parameter == "generateSamplesFile")
	generateSamplesFile = std::stoi(value);
    else if (parameter == "endEffectorFileName")
	endEffectorFileName = value;
    else if (parameter == "leftFootFileName")
	leftFootFileName = value;
    else if (parameter == "rightFootFileName")
	rightFootFileName = value;
    else if (parameter == "waistFileName")
	waistFileName = value;
    else if (parameter == "leftHipFileName")
	leftHipFileName = value;
    else if (parameter == "rightHipFileName")
	rightHipFileName = value;
    else if (parameter == "leftKneeFileName")
	leftKneeFileName = value;
    else if (parameter == "rightKneeFileName")
	rightKneeFileName = value;
    else
	std::cout << "warning: " << parameter << " is not set" << std::endl;

}
