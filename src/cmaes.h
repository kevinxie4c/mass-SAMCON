#ifndef WEIRDCMAES
#define WEIRDCMAES

#include <cmath>
#include <algorithm>
#include <random>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

class WeirdCMAES
{
public:
	Eigen::VectorXd xmean, xold;
	Eigen::VectorXd weights;
	size_t N, lambda, mu;
	double sigma, mueff;
	double cc, cs, cl, cmu, damps;
	Eigen::VectorXd pc, ps;
	Eigen::MatrixXd B, C, invsqrtC;
	Eigen::VectorXd D;
	double eigeneval, chiN;
	size_t counteval;

	std::default_random_engine gen;
	std::normal_distribution<double> dist;

	WeirdCMAES(size_t n, size_t lambda, size_t mu, double sigma, Eigen::VectorXd xmean, Eigen::VectorXd diagonal):N(n), lambda(lambda), mu(mu), sigma(sigma), D(diagonal)
	{
		using namespace Eigen;
		this->xmean = xmean;

		// selection
		weights = VectorXd(mu);
		for (size_t i = 0; i < mu; ++i)
			weights[i] = log(mu + 1.0 / 2.0) - log(i + 1); // or log((lambda + 1.0) / 2.0) - log(i + 1) ?
		weights /= weights.sum();
		mueff = pow(weights.sum(), 2) / weights.array().square().sum();
		// adaptation
		cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N);
		cs = (mueff + 2) / (N + mueff + 5);
		cl = 2 / (pow(N + 1.3, 2) + mueff);
		cmu = std::min(1 - cl, 2 * (mueff - 2 + 1 / mueff) / (pow(N + 2, 2) + mueff));
		damps = 1 + 2 * std::max(0.0, sqrt((mueff - 1) / (N + 1)) - 1) + cs;

		// dynamic
		pc = VectorXd::Zero(N);
		ps = VectorXd::Zero(N);
		B = MatrixXd::Identity(N, N);
		//D = VectorXd::Ones(N);
		C = B * D.array().square().matrix().asDiagonal() * B.transpose();
		invsqrtC = B * D.array().inverse().matrix().asDiagonal() * B.transpose();
		eigeneval = 0;
		chiN = sqrt(N) * (1.0 - 1.0 / (4.0 * N) + 1.0 / ( 21.0 * N * N));
		counteval = 0;
	}

	Eigen::VectorXd getSample()
	{
		using namespace Eigen;
		/*
		#pragma omp critical (counteval_section)
		{
			++counteval;
		}
		*/
		ArrayXd randn(N);
		for (size_t i = 0; i < N; ++i) randn[i] = dist(gen);
		return xmean + sigma * B * (D.array() * randn).matrix();
	}

	void update(Eigen::MatrixXd arx)
	{
		using namespace Eigen;
		xold = xmean;
		xmean = arx * weights;

		// assume that getSample is called sampleNum times
		counteval += lambda;

		ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * invsqrtC * (xmean - xold) / sigma;
		double hsig = ps.array().square().sum() / (1 - pow(1 - cs, 2.0 * counteval / lambda)) / N < 2 + 4 / (N + 1) ? 1 : 0;
		pc = (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mueff) * (xmean - xold) / sigma;

		MatrixXd artmp = (1.0 / sigma) * (arx - xold.replicate(1, mu));
		C = (1 - cl - cmu) * C + cl * (pc * pc.transpose() + (1 - hsig) * cc * (2 - cc) * C) + cmu * artmp * weights.asDiagonal() * artmp.transpose();

		sigma = sigma * exp((cs / damps) * (ps.norm() / chiN - 1));

		if (counteval - eigeneval > lambda / (cl + cmu) / N / 10.0)
		{
			eigeneval = counteval;
			// need enforce symmetry?
			EigenSolver<MatrixXd> es(C);
			B = es.eigenvectors().real();
			D = es.eigenvalues().real().array().sqrt();
			invsqrtC = B * D.array().inverse().matrix().asDiagonal() * B.transpose();
		}
	}
};

#endif
