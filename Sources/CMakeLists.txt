cmake_minimum_required(VERSION 3.12)
project(HIFUEvo
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


#include_directories(include)
#include_directories(external/build/include)

add_executable(cmaes ${CMAKE_SOURCE_DIR}/examples/cmaes.cpp 
                     ${CMAKE_SOURCE_DIR}/src/Optimizers/CovarianceMatrixAdaptationES.cpp
                     ${CMAKE_SOURCE_DIR}/src/FitnessFunctions/AckleyScore/AckleyScore.cpp
                     ${CMAKE_SOURCE_DIR}/src/Logging/EvoLogGenerator.cpp)
                
target_include_directories(cmaes PRIVATE ${CMAKE_SOURCE_DIR}/include
                                         ${CMAKE_SOURCE_DIR}/external/build/include)
target_link_libraries(cmaes ${CMAKE_SOURCE_DIR}/external/build/lib/libcmaespp.a
                            ${CMAKE_SOURCE_DIR}/external/build/lib/libfmt.a)