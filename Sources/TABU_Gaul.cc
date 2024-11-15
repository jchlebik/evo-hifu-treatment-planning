/**
 * @file        DoubleGaul_Tabu.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for solving real value minimization problems with GAUL library Tabu Search.
 * 
 * @version     0.5
 * 
 * @date        2019-10-20 (created) \n
 *              2020-04-04 (revised)
 */

extern "C" {
  #include <gaul.h>
}
#include <time.h>
#include "FitnessFunction.h"
#include "LoggingHooks_Gaul.h"
#include "CustomGaulHelpers.h"

/*************************************************************************************
 * GLOBALS
 * We cannot pass the run arguments through data because GAUL does not provide a way to set entity data 
 * so we have to resort to global variables
*************************************************************************************/
bool      kVerboseOutput = FALSE;     // Is the algorithm outputing data per each generation or only launch and result summaries
int       kDimension = 0;             // Dimension of the problem
double    kTopScore = INFINITY;       // If known, the optimal fitness of the problem this file is trying to solve
int       kFitnessCalls = 0;          // Number of calls to a fitness scoring function
int       kIterationsNeeded = 0;      // Number of iterations the optimizer really needed to finish
int       kPopulationId = -1;         // GAUL id of population in use
// END OF GLOBALS DECLARATION ********************************************************

/**
 * @brief calculates the overall score of a given individual.
 * 
 * @param [in, out] pop          - variable of data type representing population, defined by GAUL
 * @param [in, out] entity       - variable of data type representing currently examined individual, defined by GAUL
 * @return boolean               - should the optimization continue
 */
bool scoreFunction(population *pop,
                   entity *entity)
{
  double *chromo = (double*)(entity->chromosome[0]);
  applyConstraints(pop->len_chromosomes, &chromo);
  entity->fitness = 0.0 - fitnessFunction(pop->len_chromosomes, chromo);
  kFitnessCalls += 1;
  return true;
}
// END OF scoreFunction **************************************************************

/**
 * @brief the hook function called with the start of evaluation of every iteration.
 * 
 * @param [in]      iter          - an iteration number currently being evaluated
 * @param [in, out] best          - best entity found so far
 * @return boolean                - should the optimization continue
 */
bool iterationnHook(const int iter, 
                    entity *best)
{
  if (kVerboseOutput)
  {
    char *textToLog = createIterationLogInfo_Double(iter, best, ga_get_population_from_id(kPopulationId), kDimension);
    char *entitiesHistory = createEntitiesImprovementHistoryDump(ga_get_population_from_id(kPopulationId), best);
    if (textToLog != nullptr)
    {
      printf("%s", textToLog);
      printf("$fitnessEvaluations:%d\n", kFitnessCalls);
      //printf("%s\n", entitiesHistory);
      free(textToLog);
      free(entitiesHistory);
    }
  }
  kIterationsNeeded++;
  double *doubleCastBest = (double*)best;
  return !isConverged(iter, 
                      kFitnessCalls, 
                      kTopScore, 
                      best->fitness,
                      1,
                      &doubleCastBest);  // a pointer to a struct is a pointer to its first member so entity ** is a double **
}
// END OF iterationHook **************************************************************

/**
 * @brief Prints help for the program.
 * 
 * @return Nothing
 */
void printUsage()
{
  printf("USAGE: all arguments required\n");
  printf("./... [numberOfIterations] [dimension] [tabuListLength] [searchCount] [optimalFitness] [rngSeed] [isVerboseOutput]\n");
  printf("Constraints for optimization set in fitness function implementing source file.");
}

/**
 * @brief the main funcion of GAUL Tabu Search optimizer.
 * 
 * @param [in]      argc      - a number of arguments passed to optimizer
 * @param [in]      argv      - argumenst passed to optimizer
 * @return int                - the return value of the program
 */
int main(int argc, char **argv)
{
  population*       pop;      // Population of solutions.
  entity*           best;     // Best solution
  struct timespec   start, end;

  if (argc < 2 || argc != 8 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
  {
    printUsage();
    return EXIT_FAILURE;
  }

  int iterations = 0, listLength = 0, searchCount = 0, seed = 0;

  // cmd arguments
  iterations      = strtol(argv[1], nullptr, 10);
  kDimension      = strtol(argv[2], nullptr, 10);
  listLength      = strtol(argv[3], nullptr, 10);
  searchCount     = strtol(argv[4], nullptr, 10);
  kTopScore       = strtof(argv[5], nullptr);
  seed            = strtol(argv[6], nullptr, 10);
  kVerboseOutput  = (bool)strtol(argv[7], nullptr, 10);

  if (iterations < 1 || kDimension < 1 || listLength < 0 || searchCount < 1 )
  {
    printUsage();
    return EXIT_FAILURE;
  }

  printf("$dimensionOfProblem:%d$maximumIterations:%d$tabuListLenght:%d$searchCount:%d", 
          kDimension, iterations, listLength, searchCount);

  double *upperBounds = (double*)malloc(sizeof(double) * kDimension);
  double *lowerBounds = (double*)malloc(sizeof(double) * kDimension);
  void   *customData  = malloc(sizeof(double*) * 2);
  if (!upperBounds || !lowerBounds || ! customData)
  {
    printf("$END:errMal\n");
    return EXIT_FAILURE;
  }

  getConstraints(kDimension, &upperBounds, &lowerBounds);  //get constraints of optimized problem to preallocated arrays
  
  //set the constraints to GAUL's custom data so they can be passed to the seeding function
  ((double**)customData)[0] = lowerBounds;
  ((double**)customData)[1] = upperBounds;

  printf("$lowerBounds:");
  for (int i = 0; i < kDimension; i++)
  {
    printf("%g", lowerBounds[i]);
    if (i + 1 < kDimension)
    {
      printf(",");
    }
  }
  printf("$upperBounds:");
  for (int i = 0; i < kDimension; i++)
  {
    printf("%g", upperBounds[i]);
    if (i + 1 < kDimension)
    {
      printf(",");
    }
  }
  printf("$optimalScore:%g$rngSeed:%u$", kTopScore, seed);

  printf("\n@@@\n"); // separator of launch and run section
  random_seed(seed);

  pop = ga_genesis_double(
    0,                                  // const int              population_size
    1,                                  // const int              num_chromo
    kDimension,                         // const int              len_chromo
    nullptr,                            // GAgeneration_hook      generation_hook
    iterationnHook,                     // GAiteration_hook       iteration_hook
    nullptr,                            // GAdata_destructor      data_destructor
    nullptr,                            // GAdata_ref_incrementor data_ref_incrementor
    scoreFunction,                      // GAevaluate             evaluate
    randomUniformSeedInBounds,          // GAseed                 seed
    nullptr,                            // GAadapt                adapt
    nullptr,                            // GAselect_one           select_one
    nullptr,                            // GAselect_two           select_two
    mutateDoubleMultipointCustomStddev, // GAmutate               mutate
    nullptr,                            // GAcrossover            crossover
    nullptr,                            // GAreplace              replace
    customData                          // vpointer               User data
  );
  
  kPopulationId = ga_get_population_id(pop);

  ga_population_set_allele_mutation_prob(pop, 1.0/5.0);
  ga_population_set_allele_min_double(pop, -DBL_MAX);

  ga_population_set_tabu_parameters(
    pop,
    ga_tabu_check_double,
    listLength,
    searchCount
  );

  best = ga_get_free_entity(pop);
  ga_entity_seed(pop, best);

  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  ga_tabu(
    pop,
    best,
    iterations
  );
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);

  unsigned long long elapsedTime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
  if (kVerboseOutput)
  {
    char *entitiesHistory = createEntitiesImprovementHistoryDump(ga_get_population_from_id(kPopulationId), best);
    printf("%s\n", entitiesHistory);
  }
  const double *chromo = ((double **)best->chromosome)[0];
  printf("\n@@@\n");
  printf("$fitness:%g$time[uS]:%llu$fitnessEvaluations:%d$iterationsNeeded:%d$best:", 
          best->fitness, elapsedTime, kFitnessCalls, kIterationsNeeded - 1);
  for (int i = 0; i < kDimension; i++)
  {
    printf("%g", chromo[i]);
    if (i < kDimension - 1)
    {
      printf(",");
    }
  }
  printf("$\n");

  ga_extinction(pop);
  return EXIT_SUCCESS;
}
// END OF MAIN ***********************************************************************