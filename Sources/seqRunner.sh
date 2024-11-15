#!/bin/bash
###########################################################################################
# @file         seqRunner.sh
#
# @author       Jakub Chlebik \n
#               Faculty of Information Technology \n
#               Brno University of Technology \n
#               xchleb07@stud.fit.vutbr.com
#
# @brief        A very simple shell script to run multiple independent runs \n
#				of optimizer on a benchmarking score. Expects user to edit \n
#				the variables
#
# @version      0.1
#
# @date         2019-04-06 (created) \n
############################################################################################


ALG=$1
BENCH="Griewank"

DIM=16

for i in {1..30}
do
  digits=12       # number of digits in seed
  a=$(date +%s)
  seed=$((a*RANDOM))

  LOG_FILE="Outputs/benchmarks/${BENCH}/2000/${ALG}/${ALG}_${i}.optOut"

  while [ ${#seed} -lt 12 ]; do
      seed="${seed}$RANDOM"
  done

case $ALG in
  "ga" )
    ./ga 100 64 ${DIM} 0.8 0.1 2 0 ${seed} 1 > ${LOG_FILE}
    ;;
  "sa" )
    ./sa 2000 ${DIM} 2000 1 0 ${seed} 1 > ${LOG_FILE}
    ;;
  "de" )
    ./de 100 30 ${DIM} 1 0.8 0.5 0.7 3 1 0 ${seed} 1 > ${LOG_FILE}
     ;;
  "tabu" )
    ./tabu 1000 ${DIM} 20 5 0 ${seed} 1 > ${LOG_FILE}
     ;;
  "pso" )
    ./pso 100 75 ${DIM} 2 2 0.5 0 ${seed} 1 > ${LOG_FILE}
     ;;
  "cmaes" )
    ./cmaes 200 ${DIM} 0 ${seed} 1 > ${LOG_FILE}
     ;;
esac
done