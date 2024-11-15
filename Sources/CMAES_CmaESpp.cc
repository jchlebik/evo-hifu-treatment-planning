/**
 * @file        DoubleCMAESpp_CMAES.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for solving real value minimization problems with LibOPT Particle Swarm 
 * 
 * @version     0.5
 * 
 * @date        2019-01-09 (created) \n
 *              2019-04-04 (revised)
 */

#include <cstdlib>
#include <ctime>
#include <string>
#include <algorithm>
#include <vector>
#include <random>

#include "cmaes.h"
#include "LoggingHooks_CMAESpp.h"
#include "FitnessFunction.h"

/**************************************************************************************
 * GLOBALS
 * We cannot pass the run arguments through data so we have to resort to global variables
***************************************************************************************/
bool    kVerboseOutput = false;   // Is the algorithm outputing data per each generation or only launch and result summaries
int     kChromoLen = 0;           // Chromozome Lenght
double  kTopScore = INFINITY;     // If known, the optimal fitness of the problem this file is trying to solve
int     kFitnessCalls = 0;        // Number of calls to a fitness scoring function
int     kGenerationsNeeded = 0;   // Number of generations the optimizer really needed to finish
// END OF GLOBALS DECLARATION *********************************************************

/**
 * @brief calculates the overall score of a given individual.
 * 
 * @param [in, out] a            - agent of the swarm currently being evaluated, defined by LibOPT
 * @param [in]      arg          - list of custom arguments to be optionaly passed to fitness function, not used
 * @return doubl                 - the score given to the evaluated agent
 */
double scoreFunction(int n, 
                     double *x)
{  
  
  double *upperBound    = new double[n];
  double *lowerBound    = new double[n];
  double *mappedValues  = new double[n];

  getConstraints(n, &upperBound, &lowerBound); 

  for (int i = 0; i < n; i++)
  {
    mappedValues[i] = lowerBound[i] + (upperBound[i] - lowerBound[i]) * (1.0 - std::cos(M_PI * x[i]/10.0))/2.0;
  }

  double score = fitnessFunction(n, mappedValues);
  kFitnessCalls += 1;
  return score;
}
// END OF scoreFunction **************************************************************

/**
 * @brief the hook function called with the start of evaluation of every generation.
 * 
 * @param [in]  evo    - a class containing all cmaes data for the optimalization run
 * @param [in]  pop    - a current population
 * @return boolean     - should the optimization continue
 */
template <class T>
bool generationHook(CMAES<T> & evo)
{
  double* popFitness = (double*)evo.getPtr(CMAES<double>::FVals);
  if (kVerboseOutput)
  {
    char* textToLog = createPopulationLogInfo_Double(evo, popFitness);
    char* popDump   = createPopulationDump<double>(popFitness, evo.get(CMAES<double>::PopSize));

    if (textToLog != nullptr && popDump != nullptr)
    {
      printf("%s", textToLog);
      printf("$fitnessEvaluations:%d", kFitnessCalls);
      printf("%s\n", popDump);
      free(textToLog);
      free(popDump);
    }
    else
    {
      fprintf(stderr, "population died out.\n");
      return false;
    }
  }
  kGenerationsNeeded += 1;

  return !isConverged(evo.get(CMAES<double>::Generation), 
                     kFitnessCalls, 
                     kTopScore, 
                     evo.get(CMAES<double>::Fitness),
                     evo.get(CMAES<double>::PopSize),
                     &popFitness
  );
}
// END OF generationHook *************************************************************
/**
 * @brief runs the CMAES evolution
 * 
 * @param [in, out]    evo      - the cmaes class with all relevant information
 * @return void                 - all results of the optimization are saved and encoded in the final position of particles in the CMAES class
 */
template <class T>
void runCMAES(CMAES<T> & evo, T* arFunvals, T* upperBounds, T* lowerBounds)
{
  const double    dim = evo.get(CMAES<double>::Dimension);
  double* const*  pop;
  // Iterate until stop criterion holds
  while(!evo.testForTermination())
  {
    // Generate lambda new search points, sample population
    pop = evo.samplePopulation(); // Do not change content of pop

    /* Here you may resample each solution point pop[i] until it
       becomes feasible, e.g. for box constraints (variable
       boundaries).  
       Assumptions: the feasible domain is convex, the optimum is
       not on (or very close to) the domain boundary, initialX is
       feasible and initialStandardDeviations are sufficiently small
       to prevent quasi-infinite looping.
    */

    const double popSize = evo.get(CMAES<double>::PopSize);
    for (int i = 0; i < popSize; i++)
    {
      /*
      int geneOutOfConstraints = isInConstraints(dim, pop[i]);
      while (geneOutOfConstraints != -1)
      {
        pop = evo.reSampleGene(i, geneOutOfConstraints);
        geneOutOfConstraints = isInConstraints(dim, pop[i]);
      }
      */
      if (i < int(evo.get(CMAES<double>::Lambda)))
      {
        arFunvals[i] = scoreFunction(dim, pop[i]);
      }
    }

    // update the search distribution used for sampleDistribution()
    evo.updateDistribution(arFunvals);
    if (!generationHook<double>(evo))
    {
      break;
    }
  }
}
// END OF runCMAES *************************************************************

/**
 * @brief Prints help for the program.
 * 
 * @return Nothing.
 */
void printUsage()
{
  printf("USAGE: all arguments required\n");
  printf("./... [maxGenerations] [dimension] [optimalFitness] [rngSeed] [isVerboseOutput]\n");
  printf("starting positions and standard deviations are calculated as mean and stdev of normal distribution from the specified bounds in fitness function.\n");
}
/**
 * @brief the main funcion of GAUL Differential Evolution optimizer.
 * 
 * @param [in]      argc      - a number of arguments passed to optimizer
 * @param [in]      argv      - argumenst passed to optimizer
 * @return int                - the return value of the program
 */
int main(int argc, const char **argv)
{

  struct timespec start, end;

  if (argc < 2 || argc != 6)
  {
    printUsage();
    return EXIT_FAILURE;
  }
  
  CMAES<double> evo;
  Parameters<double> params;
  int iterations = 0, seed = 0;

  // cmd arguments
  iterations      = strtol(argv[1], nullptr, 10);
  kChromoLen      = strtol(argv[2], nullptr, 10); //dimension
  kTopScore       = strtof(argv[3], nullptr);
  seed            = strtol(argv[4], nullptr, 10);
  kVerboseOutput  = bool(strtol(argv[5], nullptr, 10));

  if (iterations < 1 || kChromoLen < 1)
  {
    printUsage();
    return EXIT_FAILURE;
  }

  double* xstart      = new double[kChromoLen];
  double* stddev      = new double[kChromoLen];
  double* upperBound  = new double[kChromoLen];
  double* lowerBound  = new double[kChromoLen];

  //srand(seed);

  getConstraints(kChromoLen, &upperBound, &lowerBound); 

  std::default_random_engine generator(seed);
 
  for(int i = 0; i < kChromoLen; i++)
  {
    std::uniform_real_distribution<double> distribution(1.0, 10.0);
    xstart[i] = distribution(generator);                                      // start at random position between bounds
    stddev[i] = (10.0 - 1.0)/3.0; 

    //std::uniform_real_distribution<double> distribution(lowerBound[i], upperBound[i]);
    //xstart[i] = rand() % int(upperBound[i] - lowerBound[i]) + lowerBound[i];  // start at random position between bounds
    //xstart[i] = distribution(generator);                                      // start at random position between bounds
    //stddev[i] = (upperBound[i] - lowerBound[i])/3.0;                          // sigma is one third of variable range
  } 

  //applyConstraints(kChromoLen, &xstart);

  params.init(kChromoLen, xstart, stddev);
  params.stopMaxIter = iterations;
  double* arFunvals = evo.init(params, seed);

  printf("$lambda:%d$genesInChromosome:%d$iterations:%d", 
          (int)evo.get(CMAES<double>::SampleSize), kChromoLen, iterations);
  printf("$lowerBounds:");
  for (int i = 0; i < kChromoLen; i++)
  {
    printf("%g", lowerBound[i]);
    if (i + 1 < kChromoLen)
    {
      printf(",");
    }
  }
  printf("$upperBounds:");
  for (int i = 0; i < kChromoLen; i++)
  {
    printf("%g", upperBound[i]);
    if (i + 1 < kChromoLen)
    {
      printf(",");
    }
  }
  printf("$optimalScore:%g$rngSeed:%u$",kTopScore, seed);
  printf("\n@@@\n"); // separator of launch and run section

  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  runCMAES<double>(evo, arFunvals, upperBound, lowerBound);
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);

  unsigned long long elapsedTime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
  printf("\n@@@\n"); // separator of run and results section
  
  double best =     evo.get(CMAES<double>::FBestEver);
  double* chromo =  evo.getNew(CMAES<double>::XBestEver); 

  printf("$fitness:%g$time[uS]:%llu$fitnessEvaluations:%d$generationsNeeded:%d$best:", 
          best, elapsedTime, kFitnessCalls, kGenerationsNeeded - 1);
  for (int i = 0; i < kChromoLen; i++)
  {
    printf("%g", chromo[i]);
    if (i < kChromoLen - 1)
    {
      printf(",");
    }
  }
  printf("$\n");
  delete[] chromo;
  delete[] xstart;
  delete[] stddev;
  delete[] upperBound;
  delete[] lowerBound;
  return EXIT_SUCCESS;
}
// END OF MAIN ***********************************************************************