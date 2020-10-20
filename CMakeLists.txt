project(mass_SAMCON)

set(CMAKE_CXX_FLAGS "-fopenmp -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
add_executable(mass.out mass.cpp BVHData.cpp MyWindow.cpp config_loader.cpp refine.cpp sample.cpp setup_frag.cpp simulator.cpp timer.cpp utility.cpp guided_samcon.cpp make_state.cpp online.cpp)
add_executable(nsr.out nsr.cpp BVHData.cpp utility.cpp config_loader.cpp)
target_link_libraries(mass.out m stdc++ dart dart-collision-ode assimp boost_system dart-gui GL GLU glut gomp pthread)
target_link_libraries(nsr.out stdc++ dart assimp boost_system)
