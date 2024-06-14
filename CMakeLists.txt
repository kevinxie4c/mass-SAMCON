cmake_minimum_required(VERSION 3.22)
project(mass_SAMCON)

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(DART REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${DART_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "-fopenmp -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
set(SRC_DIR "src")
add_executable(mass.out ${SRC_DIR}/mass.cpp ${SRC_DIR}/BVHData.cpp ${SRC_DIR}/MyWindow.cpp ${SRC_DIR}/config_loader.cpp ${SRC_DIR}/refine.cpp ${SRC_DIR}/sample.cpp ${SRC_DIR}/setup_frag.cpp ${SRC_DIR}/simulator.cpp ${SRC_DIR}/timer.cpp ${SRC_DIR}/utility.cpp ${SRC_DIR}/guided_samcon.cpp ${SRC_DIR}/make_state.cpp ${SRC_DIR}/online.cpp)
add_executable(nsr.out ${SRC_DIR}/nsr.cpp ${SRC_DIR}/BVHData.cpp ${SRC_DIR}/utility.cpp ${SRC_DIR}/config_loader.cpp)
target_link_libraries(mass.out m stdc++ dart dart-collision-ode assimp boost_system dart-gui GL GLU glut gomp pthread Eigen3::Eigen)
target_link_libraries(nsr.out stdc++ dart assimp boost_system Eigen3::Eigen)
