###########################################################################################
# @file         Makefile
#
# @author       Jakub Chlebik \n
#               Faculty of Information Technology \n
#               Brno University of Technology \n
#               xchleb07@stud.fit.vutbr.com
#
# @brief        Makefile to build optimizers
#
# @version      0.3
#
# @date         2019-10-28 (created) \n
#               2020-06-03 (revised)
############################################################################################

USER:=$(shell whoami)
CC=g++
ICC=icpc
CFLAGS=-std=c++11 -O3

LIBS_GAUL=-lgaul -lgaul_util
LIBS_LOPT=-lopt
LIBS=-lpthread -lm -ffast-math

LIBS_FFTW=-lfftw3f_omp -lfftw3f -lm -ffast-math
LIBS_HIFU_I=$(LIBS_OMP_I) -mkl -march=native -ffast-math
LIBS_OMP_I=-liomp5 -qopenmp -lpthread
LIBS_OMP_G=-liomp5 -lpthread -march=native

USR_INCPATH=/home/$(USER)/libs/include
USR_LIBPATH=/home/$(USER)/libs/lib

LIBOPTLIBFOLDERPATH=./ThirdParty/LibOPT/lib  

LOGHOOKHEADERSPATH=./LoggingFunc/
HELPERSHEADERSPATH=./Helpers/
FITNESSHEADERSPATH=./FitnessFunc/
HIFUFITNESSPATH=$(FITNESSHEADERSPATH)HIFUScore/
BENCHFITNESSPATH=$(FITNESSHEADERSPATH)Benchmarks/

CMAESPPHEADERPATH=./ThirdParty/CMA-ESpp/cma-es/
LIBOPTHEADERSPATH=./ThirdParty/LibOPT/include/

HIFUSCOREPATH=$(HIFUFITNESSPATH)heatDiffusionScore.cc
BENCHSCOREPATH=$(BENCHFITNESSPATH)GriewankScore.cc

HIFUSOURCEFILES=$(HIFUFITNESSPATH)KWaveDiffusionSolver/kWaveDiffusionSolver.cc \
								$(HIFUFITNESSPATH)utils.cc \
								$(HIFUFITNESSPATH)discreteMap.cc \
								$(HIFUFITNESSPATH)getHeatSource.cc  \
								$(HIFUFITNESSPATH)calculateHeating.cc  \
								$(HIFUFITNESSPATH)mapsFactory.cc

HELPERSPATH=Helpers/CustomGaulHelpers.cc
GAULLOGHOOKSPATH=LoggingFunc/LoggingHooks_Gaul.cc
LIBOPTLOGHOOKSPATH=LoggingFunc/LoggingHooks_LibOpt.cc
CMAESPPLOHGOOKSPATH=LoggingFunc/LoggingHooks_CMAESpp.cc

.PHONY: ga sa tabu pso cmaes ga_hifu sa_hifu tabu_hifu pso_hifu cmaes_hifu hifu_test

all: libopt opt hifu

libopt:
	@cd ./ThirdParty/LibOPT && make libopt > /dev/null

opt: 	ga de sa tabu pso cmaes

hifu: ga_hifu de_hifu sa_hifu tabu_hifu pso_hifu cmaes_hifu

hifu_test: $(HIFUSCOREPATH) $(HIFUSOURCEFILES) hifuTest.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) -I$(FITNESSHEADERSPATH) -o hifu_test

ga_hifu: $(HIFUSCOREPATH) $(HIFUSOURCEFILES) $(HELPERSPATH) $(GAULLOGHOOKSPATH) GA_Gaul.cc
		$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -I$(HELPERSHEADERSPATH) -I$(USR_INCPATH) -L$(USR_LIBPATH) -Wl,-R$(USR_LIBPATH) $(LIBS_GAUL) -o ga_hifu

de_hifu: $(HIFUSCOREPATH) $(HIFUSOURCEFILES) $(HELPERSPATH) $(GAULLOGHOOKSPATH) DE_Gaul.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -I$(HELPERSHEADERSPATH) -I$(USR_INCPATH) -L$(USR_LIBPATH) -Wl,-R$(USR_LIBPATH) $(LIBS_GAUL) -o de_hifu

sa_hifu: $(HIFUSCOREPATH) $(HIFUSOURCEFILES) $(HELPERSPATH) $(GAULLOGHOOKSPATH) SA_Gaul.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -I$(HELPERSHEADERSPATH) -I$(USR_INCPATH) -L$(USR_LIBPATH) -Wl,-R$(USR_LIBPATH) $(LIBS_GAUL) -o sa_hifu

tabu_hifu: $(HIFUSCOREPATH) $(HIFUSOURCEFILES) $(HELPERSPATH) $(GAULLOGHOOKSPATH) TABU_Gaul.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -I$(HELPERSHEADERSPATH) -I$(USR_INCPATH) -L$(USR_LIBPATH) -Wl,-R$(USR_LIBPATH) $(LIBS_GAUL) -o tabu_hifu

pso_hifu: $(HIFUSCOREPATH) $(HIFUSOURCEFILES) $(LIBOPTLOGHOOKSPATH) PSO_LibOpt.cc | libopt
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -I$(LIBOPTHEADERSPATH) -L$(LIBOPTLIBFOLDERPATH) $(LIBS_LOPT) -mkl -o pso_hifu

cmaes_hifu: $(HIFUSCOREPATH) $(HIFUSOURCEFILES) $(CMAESPPLOHGOOKSPATH) CMAES_CmaESpp.cc
	$(ICC) $(CFLAGS) $^ $(LIBS_HIFU_I) -I$(CMAESPPHEADERPATH) -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -o cmaes_hifu


ga: $(BENCHSCOREPATH) $(HELPERSPATH) $(GAULLOGHOOKSPATH) GA_Gaul.cc
	$(CC) $(CFLAGS) $^ -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -I$(HELPERSHEADERSPATH) -I$(USR_INCPATH) -L$(USR_LIBPATH) -Wl,-R$(USR_LIBPATH) $(LIBS_GAUL) $(LIBS) -o ga

de: $(BENCHSCOREPATH) $(HELPERSPATH) $(GAULLOGHOOKSPATH) DE_Gaul.cc
	$(CC) $(CFLAGS) $^ -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -I$(HELPERSHEADERSPATH) -I$(USR_INCPATH) -L$(USR_LIBPATH) -Wl,-R$(USR_LIBPATH) $(LIBS_GAUL) $(LIBS) -o de

sa: $(BENCHSCOREPATH) $(HELPERSPATH) $(GAULLOGHOOKSPATH) SA_Gaul.cc
	$(CC) $(CFLAGS) $^ -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -I$(HELPERSHEADERSPATH) -I$(USR_INCPATH) -L$(USR_LIBPATH) -Wl,-R$(USR_LIBPATH) $(LIBS_GAUL) $(LIBS) -o sa

tabu: $(BENCHSCOREPATH) $(HELPERSPATH) $(GAULLOGHOOKSPATH) TABU_Gaul.cc
	$(CC) $(CFLAGS) $^ -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -I$(HELPERSHEADERSPATH) -I$(USR_INCPATH) -L$(USR_LIBPATH) -Wl,-R$(USR_LIBPATH) $(LIBS_GAUL) $(LIBS) -o tabu

pso: $(BENCHSCOREPATH) $(LIBOPTLOGHOOKSPATH) PSO_LibOpt.cc | libopt
	$(ICC) $(CFLAGS) $^ -I$(LIBOPTHEADERSPATH) -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -L $(LIBOPTLIBFOLDERPATH) $(LIBS_LOPT) -o pso

cmaes: $(BENCHSCOREPATH) $(CMAESPPLOHGOOKSPATH) CMAES_CmaESpp.cc
	$(CC) $(CFLAGS) $^ -I$(CMAESPPHEADERPATH) -I$(FITNESSHEADERSPATH) -I$(LOGHOOKHEADERSPATH) -lm -o cmaes

clean:
	rm ga | rm sa | rm tabu | rm pso | rm de | rm cmaes | rm ga_hifu | rm sa_hifu | rm tabu_hifu | rm pso_hifu | rm de_hifu | rm cmaes_hifu | rm hifu_test
