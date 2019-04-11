darttest.out: darttest.cpp BVHData.o BVHData.h cmaes.h
	g++ -g -O3 -std=c++11 -fopenmp -o darttest.out darttest.cpp BVHData.o -ldart -ldart-collision-ode -lassimp -lboost_system -ldart-gui -lGL -lGLU -lglut
BVHData.o: BVHData.cpp BVHData.h
	g++ -g -O3 -c -std=c++11 -o BVHData.o BVHData.cpp
