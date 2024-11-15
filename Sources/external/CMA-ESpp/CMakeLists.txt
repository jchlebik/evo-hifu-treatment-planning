cmake_minimum_required(VERSION 3.5)
project(cmaespp
        LANGUAGES CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

  set(CMAKE_CXX_EXTENSIONS OFF)
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    #include_directories(cma-es)

# doxygen, source:
# http://permalink.gmane.org/gmane.comp.programming.tools.cmake.user/4930
find_package(Doxygen)
  if(DOXYGEN_FOUND)
  set(DOXYGEN_INPUT doxygen.cfg)
  add_custom_target(apidoc
    COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
    COMMAND ${CMAKE_COMMAND} -E echo "Done."
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endif(DOXYGEN_FOUND)



#add_executable(evo2 example2.cpp)
endif()

add_subdirectory(cma-es)

set(HEADER_LIST "cma-es/include/cmaes.h" "cma-es/include/parameters.h" "cma-es/include/random.h" "cma-es/include/timings.h" "cma-es/include/utils.h")

install(TARGETS cmaespp DESTINATION lib)
install(FILES ${HEADER_LIST} DESTINATION include/cmaespp)
#add_executable(evo1 example1.cpp)
#target_compile_features(evo1 PUBLIC cxx_std_17)
#target_link_libraries(evo1 PUBLIC cmaespp)
#target_include_directories(evo1 PUBLIC cma-es)
