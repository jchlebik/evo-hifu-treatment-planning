/**
 * @file        hifuTest.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       A small program capable of running one heat diffusion simulation \n
 *              and printing result to the standard error output
 * 
 * @version     0.1
 * 
 * @date        2020-04-06 (created) \n
 */

#include <iostream>
#include "FitnessFunc/FitnessFunction.h"
#include <ctime>

/**
 * @brief the main funcion of small HIFU simulation test.
 * 
 * @return int    - the return value of the program
 */
int main()
{

  struct timespec   start, end;
  double arr[] = {219.612,276.122,9.565,89.7943,217.881,219.77,8.39193,87.1744,228.734,256.377,7.65081,83.1919,256.724,276.753,9.24244,96.2361,239.166,281.089,7.07247,77.2956,218.265,216.67,7.45452,48.7604,217.265,248.73,8.93117,59.4931,271.995,255.494,9.14906,91.7417,238.4,208.947,7.02642,73.5941,269.383,214.727,9.09386,60.7772};

  int n = 40;

	/*
  double* upperBound = (double*)malloc(n*sizeof(double));
  double* lowerBound = (double*)malloc(n*sizeof(double));
  getConstraints(n, &upperBound, &lowerBound); 
  applyConstraints(n, &arr);
	*/


  clock_gettime(CLOCK_MONOTONIC_RAW, &start);

  int res = fitnessFunction(n, arr);
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);
  unsigned long long elapsedTime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
  std::cerr << "Run time: " << elapsedTime << "[uS]\n";
  std::cerr << "Resulting fitness: " << res << "\n";
  return 0;
}
