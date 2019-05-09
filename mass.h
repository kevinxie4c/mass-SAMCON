#ifndef GUIDED_H
#define GUIDED_H

#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Core>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#endif
#include "control_fragment.h"
#include "simulator.h"

extern size_t numFrag;
extern std::string taskFileName;
extern std::vector<size_t> walk;
extern std::vector<ControlFragment> frags;
extern std::vector<Eigen::VectorXd> initMean;
extern std::vector<Simulator> simulators;

void setUpFrags(bool useMass);

void refine(bool useMass);

#endif
