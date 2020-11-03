#ifndef CONFIG_LOADER
#define CONFIG_LOADER

#include <string>

namespace Config
{
	// Simulation
	extern double frictionCoeff;
	extern double restitutionCoeff;
	extern double gravity;
	extern double ERP;
	extern double CFM;
	extern double timeStep;
	extern double frameTime;
	extern size_t stepPerFrame;

	// cost function
	extern double wp, wr, we, wb, w_zmp;

	// file name
	extern std::string bvhFileName;
	extern std::string geometryConfigFileName;
	extern std::string hingeJointListFileName;
	extern std::string freeNodesListFileName;
	extern std::string outputFileName;
	extern std::string basedir;
	extern std::string scalesFileName;
	extern std::string covFileName;
	extern std::string stiffnessFileName;
	extern std::string dampingFileName;
	extern std::string massFileName;
	extern std::string initMeanFileName;
	extern std::string endEffectorFileName;
	extern std::string leftFootFileName;
	extern std::string rightFootFileName;
	extern std::string waistFileName;
	extern std::string leftHipFileName;
	extern std::string rightHipFileName;
	extern std::string leftKneeFileName;
	extern std::string rightKneeFileName;
	extern std::string initStateFileName;
	extern std::string extForceFileName;

	// algorithm
	extern double scale;
	extern size_t startFrame;
	extern size_t endFrame;
	extern double groundOffset;
	extern size_t sampleNum, saveNum;
	extern double percentageDiscard;
	extern int convexExp;
	extern size_t groupNum;
	extern size_t slidingWindow;
	extern size_t updateWindow;
	extern size_t trialMin, trialMax;
	extern size_t notImproveMax;
	extern double goodEnough, failThreshold;
	extern bool stablePD;
	extern bool zeroInitialVelocities;
	extern bool zeroInitialAngularVelocities;
	extern bool useSampleNumAsLambda;
	extern bool useAFT; // we don't use it now
	extern double initSigma;
	extern size_t maxFailure;
	extern size_t rank;
	extern size_t dRank;
	extern double scaleMassMatrix;
	extern bool useEigenvalueScale;
	extern bool autoScaleMassMatrix;
	extern bool useInverseForce;
	extern bool useEigenvector;
	extern bool useEigenvalue;
	extern std::string flexibleJointsFileName;
	extern std::string forceFileName;
	extern std::string noIDJointsFileName;
	extern bool useCompensator;
	extern bool k_cmp;
	extern bool d_cmp;
	extern double inverseForceRatio;
	extern bool originalSAMCON;

	// others
	extern bool showWindow;
	extern bool onlyLogAndFinal;
	extern size_t loopNum;
	extern bool generateSamplesFile;


	void load(const std::string &filename);

	void setParameter(const std::string &parameter, const std::string &value);

};

#endif
