darttest.out: darttest.cpp BVHData.o BVHData.h cmaes.h
	g++ -std=c++11 -fopenmp -o darttest.out darttest.cpp BVHData.o -ldart -ldart-collision-ode
BVHData.o: BVHData.cpp BVHData.h
	g++ -c -std=c++11 -o BVHData.o BVHData.cpp
