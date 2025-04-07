cmake_minimum_required(VERSION 3.10)
project(SimplexSolver)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)





add_executable(SimplexSolver 
    src/main.cpp 
    src/matrix.cpp 
    src/SimplexEngine.cpp 
    src/SimplexSolver.cpp
)

target_include_directories(SimplexSolver PRIVATE headers)
