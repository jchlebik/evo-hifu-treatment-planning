/**
 * @file        DoubleLibOpt_PSO.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for solving real value minimization problems with LibOPT Particle Swarm 
 * 
 * @version     0.6
 * 
 * @date        2019-01-09 (created) \n
 *              2020-04-04 (revised)
 */

#include <cstdlib>
#include <ctime>
#include <string>
#include <algorithm>
#include <vector>
#include "FitnessFunction.h"
#include "LoggingHooks_LibOpt.h"
extern "C"{
  #include "LibOpt_PSO.h"
}

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
double scoreFunction(Agent *a, 
                     va_list arg)
{
  double score = fitnessFunction(a->n, a->x);
  kFitnessCalls += 1;
  return score;
}
// END OF scoreFunction **************************************************************

/**
 * @brief the hook function called with the start of evaluation of every generation.
 * 
 * @param [in]      generation    - a generation number currently being evaluated
 * @param [in, out] searchSpace   - a structure containing the data of the search space and agents in it
 * @return boolean                - should the optimization continue
 */
bool generationHook(const int generation, SearchSpace *searchSpace)
{

  if (kVerboseOutput)
  {
    char* textToLog = createPopulationLogInfo_Double(generation, searchSpace);
    char* popDump = createPopulationDump(searchSpace);

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

  return !isConverged(generation, 
                      kFitnessCalls,
                      kTopScore,
                      searchSpace->a[searchSpace->best]->fit,
                      searchSpace->m,
                      (double**)searchSpace->a
  );
}
// END OF generationHook *************************************************************

/**
 * @brief Prints help for the program.
 * 
 * @return Nothing.
 */
void printUsage()
{
  printf("USAGE: all arguments required\n");
  printf("./... [generationsNum] [populationSize] [genesInChromosome] [c1] [c2] [w] [optimalFitness] [rngSeed] [isVerboseOutput]\n");
  printf("Constraints for optimization set in fitness function implementing source file.\n");
}

/**
 * @brief function to create and fill LibOPT SearchSpace structure
 * 
 * @param [in]      argv      - arguments passed to the optimizer
 * @return SearchSpace        - structure containing all relevant information about the SearchSpace needed to run PSO
 */
SearchSpace* prepareSearchSpace(const char** argv)
{
  SearchSpace *s; 

  int generations = 0, popSize = 0, seed = 0;
  float c1 = 0.0, c2 = 0.0, w = 0.0;

  // cmd arguments
  generations     = strtol(argv[1], nullptr, 10);
  popSize         = strtol(argv[2], nullptr, 10);
  kChromoLen      = strtol(argv[3], nullptr, 10);
  c1              = strtof(argv[4], nullptr);
  c2              = strtof(argv[5], nullptr);
  w               = strtof(argv[6], nullptr);
  //w_min           = strtof(argv[7], nullptr);
  //w_max           = strtof(argv[8], nullptr);
  kTopScore       = strtof(argv[7], nullptr);
  seed            = strtol(argv[8], nullptr, 10);
  kVerboseOutput  = (bool)strtol(argv[9], nullptr, 10);

  if (generations < 1 || popSize < 1 || kChromoLen < 1 || c1 < 0 || c2 < 0 || w < 0)
  {
    printUsage();
    return nullptr;
  }

  printf("$populationSize:%d$genesInChromosome:%d$generationsNum:%d", popSize, kChromoLen, generations);
  printf("$c1:%g$c2:%g$w:%g", c1, c2, w);

  s = CreateSearchSpace(popSize, kChromoLen, _PSO_);
  s->iterations = generations;
  s->c1 = c1;
  s->c2 = c2;
  s->w = w;
  //s->w_min = w_min;
  //s->w_max = w_max;

  getConstraints(s->n, &(s->UB), &(s->LB));

  printf("$lowerBounds:");
  for (int i = 0; i < kChromoLen; i++)
  {
    printf("%g", s->LB[i]);
    if (i + 1 < kChromoLen)
    {
      printf(",");
    }
  }
  printf("$upperBounds:");
  for (int i = 0; i < kChromoLen; i++)
  {
    printf("%g", s->UB[i]);
    if (i + 1 < kChromoLen)
    {
      printf(",");
    }
  }
  printf("$optimalScore:%g$rngSeed:%u$", kTopScore, seed);
  printf("\n@@@\n"); // separator of launch and run section
  srandinter(seed);
  return s;
}
// END OF prepareSearchSpace *************************************************************

/**
 * @brief runs the PSO with Adaptive innertia weights
 * 
 * @param [in, out]    s      - a search space in which the particles exists
 * @return void               - all results of the optimization are saved and encoded in the final position of particles in search space
 */
void runPSO_Custom(SearchSpace *s, ...)
{
  double prob;
  va_list arg, argtmp;
  va_start(arg, s);
  va_copy(argtmp, arg);
  //printf("called pso run\n");

  for (int i = 0; i < s->m; i++)
  {
    //printf("constraint application %d\n", i);
    applyConstraints(s->n, &(s->a[i]->x));
  }

  EvaluateSearchSpace(s, _PSO_, scoreFunction, arg); /* Initial evaluation */
  // Adaptive innertia weights
  //for (int i = 0; i < s->m; i++)
  //{
  //  s->a[i]->pfit = s->a[i]->fit;
  //}
  if (!generationHook(0, s)) 
  {
    va_end(arg);
    return;
  }
  for (int t = 1; t <= s->iterations; t++)
  {
    va_copy(arg, argtmp);
    /* for each particle */
    for (int i = 0; i < s->m; i++)
    {
      UpdateParticleVelocity(s, i);
      UpdateParticlePosition(s, i);
      applyConstraints(s->n, &(s->a[i]->x));
      //CheckAgentLimits(s, s->a[i]);
    }    
    
    EvaluateSearchSpace(s, _PSO_, scoreFunction, arg);

    //Adaptive innertia weights part
    //prob = ComputeSuccess(s);                       /* Equations 17 and 18 */
    //s->w = (s->w_max - s->w_min) * prob + s->w_min; /* Equation 20 */
    //for (int i = 0; i < s->m; i++)
    //{
    //  s->a[i]->pfit = s->a[i]->fit;
    //}

    va_copy(arg, argtmp); 
    if (!generationHook(t, s)) 
    {
      break;
    }
  }
  va_end(arg);
}
// END OF runPSO_Custom *************************************************************

/**
 * @brief the main funcion of GAUL Differential Evolution optimizer.
 * 
 * @param [in]      argc      - a number of arguments passed to optimizer
 * @param [in]      argv      - argumenst passed to optimizer
 * @return int                - the return value of the program
 */
int main(int argc, const char **argv)
{
  SearchSpace *s = nullptr;
  prtFun Evaluate;

  struct timespec start, end;

  if (argc < 2 || argc != 10 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
  {
    printUsage();
    return EXIT_FAILURE;
  }

  s = prepareSearchSpace(argv);
  if (!s)
  {
    return EXIT_FAILURE;
  }

  InitializeSearchSpace(s, _PSO_);

  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  runPSO_Custom(s);
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);

  unsigned long long elapsedTime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

  printf("\n@@@\n"); // separator of run and results section
  Agent* best = s->a[s->best];
  const double *chromo = (best->x);
  printf("$fitness:%g$time[uS]:%llu$fitnessEvaluations:%d$generationsNeeded:%d$best:", 
          best->fit, elapsedTime, kFitnessCalls, kGenerationsNeeded - 1);
  for (int i = 0; i < kChromoLen; i++)
  {
    printf("%g", chromo[i]);
    if (i < kChromoLen - 1)
    {
      printf(",");
    }
  }
  printf("$\n");
  DestroySearchSpace(&s, _PSO_);
  
  return EXIT_SUCCESS;
}
// END OF MAIN ***********************************************************************