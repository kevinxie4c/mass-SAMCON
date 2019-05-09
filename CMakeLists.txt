project(mass_SAMCON)

file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
set (CMAKE_CXX_FLAGS "-fopenmp -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
add_executable(mass.out ${SRC_FILES})
target_link_libraries(mass.out m stdc++ dart dart-collision-ode assimp boost_system dart-gui GL GLU glut gomp pthread)
