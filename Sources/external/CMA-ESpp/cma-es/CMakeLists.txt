
add_library(cmaespp STATIC cmaes.cpp parameters.cpp random.cpp timings.cpp include/utils.h)
target_include_directories(cmaespp PUBLIC include)

target_compile_features(cmaespp PRIVATE cxx_std_17)