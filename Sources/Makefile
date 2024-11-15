################################################################################
# @file         Makefile                                                       #
#                                                                              #
# @author       Jakub Chlebik \n                                               #
#               Faculty of Information Technology \n                           #
#               Brno University of Technology \n                               #
#               xchleb07@stud.fit.vutbr.com                                    #
#                                                                              #
# @brief        Makefile to build optimizers                                   #
#                                                                              #
# @version      0.3                                                            #
#                                                                              #
# @date         2019-10-28 (created) \n                                        #
#               2020-06-03 (revised)                                           #
################################################################################

USER:=$(shell whoami)
CC=g++
ICC=icpc
CFLAGS=-std=c++11 -O3

# Necessary libraries
LIBS_GAUL=-lgaul -lgaul_util
LIBS_LOPT=-lopt
LIBS=-lpthread -lm -ffast-math

LIBS_FFTW=-lfftw3f_omp -lfftw3f -lm -ffast-math
LIBS_HIFU_I=$(LIBS_OMP_I) -mkl -march=native -ffast-math
LIBS_OMP_I=-liomp5 -qopenmp -lpthread
LIBS_OMP_G=-liomp5 -lpthread -march=native


#Please, Fill in GAUL paths
GAUL_INC_PATH=./ThirdParty/gaul-devel-0.1850-0/Build/include/
GAUL_LIB_PATH=./ThirdParty/gaul-devel-0.1850-0/Build/lib/

# Internal paths
LIB_OPT_LIB_PATH=./ThirdParty/LibOPT/lib/
LOG_HOOK_INC_PATH=./LoggingFunc/
LOG_HELPERS_INC_PATH=./Helpers/
FITNESS_INC_PATH=./FitnessFunc/

# Fitness functions paths
HIFU_FITNESS_PATH=$(FITNESS_INC_PATH)HIFUScore/
BENCH_FITNESS_PATH=$(FITNESS_INC_PATH)Benchmarks/

# Optimization algorithm paths
CMAESPP_INC_PATH=./ThirdParty/CMA-ESpp/cma-es/
LIB_OPT_INC_PATH=./ThirdParty/LibOPT/include/

HIFU_SCORE_PATH=$(HIFU_FITNESS_PATH)heatDiffusionScore.cc
BENCH_SCORE_PATH=$(BENCH_FITNESS_PATH)GriewankScore.cc

HIFU_SOURCE_FILES=$(HIFU_FITNESS_PATH)KWaveDiffusionSolver/kWaveDiffusionSolver.cc \
                  $(HIFU_FITNESS_PATH)utils.cc \
                  $(HIFU_FITNESS_PATH)discreteMap.cc \
                  $(HIFU_FITNESS_PATH)getHeatSource.cc  \
                  $(HIFU_FITNESS_PATH)calculateHeating.cc  \
                  $(HIFU_FITNESS_PATH)mapsFactory.cc

HELPER_SOURCE=Helpers/CustomGaulHelpers.cc
GAUL_LOG_HOOK_SOURCE=LoggingFunc/LoggingHooks_Gaul.cc
LIB_OPT_LOG_HOOK_SOURCE=LoggingFunc/LoggingHooks_LibOpt.cc
CMAESPP_LOG_HOOK_SOURCE=LoggingFunc/LoggingHooks_CMAESpp.cc


################################################################################
#                                 Build                                        #
################################################################################

.PHONY: ga sa tabu pso cmaes ga_hifu sa_hifu tabu_hifu pso_hifu cmaes_hifu hifu_test

all: libopt opt hifu

libopt:
	@cd ./ThirdParty/LibOPT && make libopt > /dev/null

opt: 	ga de sa tabu pso cmaes

hifu: ga_hifu de_hifu sa_hifu tabu_hifu pso_hifu cmaes_hifu

hifu_test: $(HIFU_SCORE_PATH) $(HIFU_SOURCE_FILES) hifuTest.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) -I$(FITNESS_INC_PATH) \
	-o hifu_test

ga_hifu: $(HIFU_SCORE_PATH) $(HIFU_SOURCE_FILES) $(HELPER_SOURCE) $(GAUL_LOG_HOOK_SOURCE) GA_Gaul.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) \
	-I$(FITNESS_INC_PATH) -I$(LOG_HOOK_INC_PATH) \
	-I$(LOG_HELPERS_INC_PATH) -I$(GAUL_INC_PATH) \
	-L$(GAUL_LIB_PATH) -Wl,-R$(GAUL_LIB_PATH) $(LIBS_GAUL) \
	-o ga_hifu

de_hifu: $(HIFU_SCORE_PATH) $(HIFU_SOURCE_FILES) $(HELPER_SOURCE) $(GAUL_LOG_HOOK_SOURCE) DE_Gaul.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) \
	-I$(FITNESS_INC_PATH) -I$(LOG_HOOK_INC_PATH) -I$(LOG_HELPERS_INC_PATH) \
	-I$(GAUL_INC_PATH) \
	-L$(GAUL_LIB_PATH) -Wl,-R$(GAUL_LIB_PATH) $(LIBS_GAUL) \
	-o de_hifu

sa_hifu: $(HIFU_SCORE_PATH) $(HIFU_SOURCE_FILES) $(HELPER_SOURCE) $(GAUL_LOG_HOOK_SOURCE) SA_Gaul.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) \
	-I$(FITNESS_INC_PATH) -I$(LOG_HOOK_INC_PATH) -I$(LOG_HELPERS_INC_PATH) \
	-I$(GAUL_INC_PATH) \
	-L$(GAUL_LIB_PATH) -Wl,-R$(GAUL_LIB_PATH) $(LIBS_GAUL) \
	-o sa_hifu

tabu_hifu: $(HIFU_SCORE_PATH) $(HIFU_SOURCE_FILES) $(HELPER_SOURCE) $(GAUL_LOG_HOOK_SOURCE) TABU_Gaul.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) \
	-I$(FITNESS_INC_PATH) -I$(LOG_HOOK_INC_PATH) -I$(LOG_HELPERS_INC_PATH) \
	-I$(GAUL_INC_PATH) -L$(GAUL_LIB_PATH) \
	-Wl,-R$(GAUL_LIB_PATH) $(LIBS_GAUL) \
	-o tabu_hifu

pso_hifu: $(HIFU_SCORE_PATH) $(HIFU_SOURCE_FILES) $(LIB_OPT_LOG_HOOK_SOURCE) PSO_LibOpt.cc | libopt
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) \
	-I$(FITNESS_INC_PATH) -I$(LOG_HOOK_INC_PATH) -I$(LIB_OPT_INC_PATH) \
	-L$(LIB_OPT_LIB_PATH) $(LIBS_LOPT) \
	-mkl -o pso_hifu

cmaes_hifu: $(HIFU_SCORE_PATH) $(HIFU_SOURCE_FILES) $(CMAESPP_LOG_HOOK_SOURCE) CMAES_CmaESpp.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) -I$(CMAESPP_INC_PATH) \
	-I$(FITNESS_INC_PATH) -I$(LOG_HOOK_INC_PATH) \
	-o cmaes_hifu


ga: $(BENCH_SCORE_PATH) $(HELPER_SOURCE) $(GAUL_LOG_HOOK_SOURCE) GA_Gaul.cc
	$(CC) $(CFLAGS) $^ -I$(FITNESS_INC_PATH) -I$(LOG_HOOK_INC_PATH) \
	-I$(LOG_HELPERS_INC_PATH) -I$(GAUL_INC_PATH) \
	-L$(GAUL_LIB_PATH) -Wl,-R$(GAUL_LIB_PATH) $(LIBS_GAUL) $(LIBS) \
	-o ga

de: $(BENCH_SCORE_PATH) $(HELPER_SOURCE) $(GAUL_LOG_HOOK_SOURCE) DE_Gaul.cc
	$(CC) $(CFLAGS) $^ -I$(FITNESS_INC_PATH) -I$(LOG_HOOK_INC_PATH) \
	-I$(LOG_HELPERS_INC_PATH) -I$(GAUL_INC_PATH) \
	-L$(GAUL_LIB_PATH) -Wl,-R$(GAUL_LIB_PATH) $(LIBS_GAUL) $(LIBS) \
	-o de

sa: $(BENCH_SCORE_PATH) $(HELPER_SOURCE) $(GAUL_LOG_HOOK_SOURCE) SA_Gaul.cc
	$(CC) $(CFLAGS) $^ -I$(FITNESS_INC_PATH) -I$(LOG_HOOK_INC_PATH) \
	-I$(LOG_HELPERS_INC_PATH) -I$(GAUL_INC_PATH) \
	-L$(GAUL_LIB_PATH) -Wl,-R$(GAUL_LIB_PATH) $(LIBS_GAUL) $(LIBS) \
	-o sa

tabu: $(BENCH_SCORE_PATH) $(HELPER_SOURCE) $(GAUL_LOG_HOOK_SOURCE) TABU_Gaul.cc
	$(CC) $(CFLAGS) $^ -I$(FITNESS_INC_PATH) -I$(LOG_HOOK_INC_PATH) \
	-I$(LOG_HELPERS_INC_PATH) -I$(GAUL_INC_PATH) -L$(GAUL_LIB_PATH) \
	-Wl,-R$(GAUL_LIB_PATH) $(LIBS_GAUL) $(LIBS) \
	-o tabu

pso: $(BENCH_SCORE_PATH) $(LIB_OPT_LOG_HOOK_SOURCE) PSO_LibOpt.cc | libopt
	$(ICC) $(CFLAGS) $^ -I$(LIB_OPT_INC_PATH) -I$(LOG_HOOK_INC_PATH) \
	-I$(FITNESS_INC_PATH) -I$(LOG_HELPERS_INC_PATH) \
	-L $(LIB_OPT_LIB_PATH) $(LIBS_LOPT) \
	-o pso

cmaes: $(BENCH_SCORE_PATH) $(CMAESPP_LOG_HOOK_SOURCE) CMAES_CmaESpp.cc
	$(CC) $(CFLAGS) $^ -I$(CMAESPP_INC_PATH) -I$(LOG_HOOK_INC_PATH) \
	-I$(FITNESS_INC_PATH) \
	-I$(LOG_HELPERS_INC_PATH) -lm \
	-o cmaes

clean:
	rm -f ga 
	rm -f sa
	rm -f tabu
	rm -f pso
	rm -f de
	rm -f cmaes
	rm -f ga_hifu
	rm -f sa_hifu
	rm -f tabu_hifu
	rm -f pso_hifu
	rm -f de_hifu
	rm -f cmaes_hifu
	rm -f hifu_test
