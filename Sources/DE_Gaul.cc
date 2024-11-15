/**
 * @file        DoubleGaul_DE.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for solving real value minimization problems with GAUL library Differential Evolution
 * 
 * @version     0.4
 * 
 * @date        2019-01-05 (created) \n
 *              2020-04-04 (revised)
 */

extern "C" {
  #include <gaul.h>
}
#include <time.h>
#include "FitnessFunction.h"
#include "LoggingHooks_Gaul.h"
#include "CustomGaulHelpers.h"

/**************************************************************************************
 * GLOBALS
 * We cannot pass the run arguments through data because GAUL does not provide a way to set entity data 
 * so we have to resort to global variables
***************************************************************************************/
bool    kVerboseOutput = FALSE;   // Is the algorithm outputing data per each generation or only launch and result summaries
int     kChromoLen = 0;           // Chromozome Lenght
double  kTopScore = INFINITY;     // If known, the optimal fitness of the problem this file is trying to solve
int     kFitnessCalls = 0;        // Number of calls to a fitness scoring function
int     kGenerationsNeeded = 0;   // Number of generations the optimizer really needed to finish
// END OF GLOBALS DECLARATION *********************************************************

/**
 * @brief calculates the overall score of a given individual.
 * 
 * @param [in]      pop          - variable of data type representing population, defined by GAUL
 * @param [in, out] entity       - variable of data type representing currently examined individual, defined by GAUL
 * @return boolean               - should the optimization continue
 */
bool scoreFunction(population *pop, 
                   entity *entity)
{
  double *chromo = (double*)(entity->chromosome[0]);
  applyConstraints(pop->len_chromosomes, &chromo);
  double score = fitnessFunction(pop->len_chromosomes, chromo);
  entity->fitness = -score;
  kFitnessCalls += 1;
  return true;
}
// END OF scoreFunction **************************************************************

/**
 * @brief the hook function called with the start of evaluation of every generation.
 * 
 * @param [in]      generation    - a generation number currently being evaluated
 * @param [in, out] pop           - variable of data type representing population, defined by GAUL
 * @return boolean                - should the optimization continue
 */
bool generationHook(const int generation, 
                    population *pop)
{
  if (kVerboseOutput)
  {
    char *textToLog = createPopulationLogInfo_Double(generation, pop);
    char *popDump = createPopulationDump(pop);
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
  entity *best = ga_get_entity_from_rank(pop, 0);

  kGenerationsNeeded += 1;
  // Stopping criterion, hook returning false stops the optimization process
  return !isConverged(generation, 
                      kFitnessCalls, 
                      kTopScore, 
                      best->fitness,
                      pop->size,
                      (double**)pop->entity_iarray);  // a pointer to a struct is a pointer to its first member so entity ** is a double **
}
// END OF generationHook *************************************************************

/**
 * @brief Prints help for the program.
 * 
 * @return Nothing
 */
void printUsage()
{
  const int strategyTypeIds[] = {0, 1, 2, 3};
  const int crossoverTypeIds[] = {0, 1, 2};
  printf("USAGE: all arguments required\n");
  printf("./... [generationsNum] [populationSize] [genesInChromosome] [numPerturbed] [crossoverFactor] [weightMin] [weightMax] [strategyType] [crossoverType] [optimalFitness] [rngSeed] [isVerboseOutput]\n");
  printf("Constraints for optimization set in fitness function implementing source file.\n");
  printf("Strategy types:\n");
  for (int i = 1; i < 4; i++)
  {
    printf("\t%d - %s", strategyTypeIds[i], getStrategyTypeByName((ga_de_strategy_type)strategyTypeIds[i]));
    if (strategyTypeIds[i] == 0)
    {
      printf(" <----- DOES NOT WORK");
    }
    printf("\n");
  }
  printf("Crossover types:\n");
  for (int i = 1; i < 3; i++)
  {
    printf("\t%d - %s", crossoverTypeIds[i], getCrossoverTypeByName((ga_de_crossover_type)crossoverTypeIds[i]));
    if (crossoverTypeIds[i] == 0)
    {
      printf(" <----- DOES NOT WORK");
    }
    printf("\n");
  }
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
  population *pop = nullptr;    // Population of solutions.
  const entity *best = nullptr; // Best entity found after GA.
  struct timespec start, end;

  if (argc < 2 || argc != 13 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
  {
    printUsage();
    return EXIT_FAILURE;
  }
  int generations = 0, popSize = 0, numPerturbed = 0, seed = 0;
  float crossoverFactor = 0.0, weightMin = 0.0, weightMax = 0.0;
  ga_de_strategy_type strategyType;
  ga_de_crossover_type crossoverType;

  // cmd arguments
  generations       = strtol(argv[1], nullptr, 10);
  popSize           = strtol(argv[2], nullptr, 10);
  kChromoLen        = strtol(argv[3], nullptr, 10);
  numPerturbed      = strtol(argv[4], nullptr, 10);
  crossoverFactor   = strtof(argv[5], nullptr);
  weightMin         = strtof(argv[6], nullptr);
  weightMax         = strtof(argv[7], nullptr);
  strategyType      = (ga_de_strategy_type)strtol(argv[8], nullptr, 10);
  crossoverType     = (ga_de_crossover_type)strtol(argv[9], nullptr, 10);
  kTopScore         = strtol(argv[10], nullptr, 10);
  seed              = strtol(argv[11], nullptr, 10);
  kVerboseOutput    = (bool)strtol(argv[12], nullptr, 10);

  if (strategyType < 1 || strategyType > 3 || 
      crossoverType < 1 || crossoverType > 2 ||
      generations < 1 || popSize < 4 ||
      kChromoLen < 1 || numPerturbed < 1 || 
      crossoverFactor < 0 || crossoverFactor > 1 ||
      weightMin < 0 || weightMax < 0 ||
      weightMax < weightMin)
  {
    printUsage();
    return EXIT_FAILURE;
  }

  printf("$populationSize:%d$genesInChromosome:%d$generationsNum:%d", 
          popSize, kChromoLen, generations);
  printf("$numPerturbed:%d$crossoverFactor:%g$weightMin:%g$weightMax:%g$strategyType:%d$crossoverType:%d$strategyTypeName:%s$crossoverTypeName:%s", 
          numPerturbed, crossoverFactor, weightMin, weightMax, strategyType, crossoverType, 
          getStrategyTypeByName(strategyType), getCrossoverTypeByName(crossoverType));
  
  double *upperBounds = (double*)malloc(sizeof(double) * kChromoLen);
  double *lowerBounds = (double*)malloc(sizeof(double) * kChromoLen);
  void   *customData  = malloc(sizeof(double*) * 2);
  if (!upperBounds || !lowerBounds || ! customData)
  {
    printf("$END:errMal\n");
    return EXIT_FAILURE;
  }

  getConstraints(kChromoLen, &upperBounds, &lowerBounds);  //get constraints of optimized problem to preallocated arrays
  
  //set the constraints to GAUL's custom data so they can be passed to the seeding function
  ((double**)customData)[0] = lowerBounds;
  ((double**)customData)[1] = upperBounds;

  printf("$lowerBounds:");
  for (int i = 0; i < kChromoLen; i++)
  {
    printf("%g", lowerBounds[i]);
    if (i + 1 < kChromoLen)
    {
      printf(",");
    }
  }
  printf("$upperBounds:");
  for (int i = 0; i < kChromoLen; i++)
  {
    printf("%g", upperBounds[i]);
    if (i + 1 < kChromoLen)
    {
      printf(",");
    }
  }
  printf("$optimalScore:%g$rngSeed:%u$", kTopScore, seed);
  printf("\n@@@\n"); // separator of launch and run section
  random_seed(seed);

  pop = ga_genesis_double(
    popSize,                    // const int              population_size
    1,                          // const int              num_chromo
    kChromoLen,                 // const int              len_chromo
    generationHook,             // GAgeneration_hook      generation_hook
    nullptr,                    // GAiteration_hook       iteration_hook
    nullptr,                    // GAdata_destructor      data_destructor
    nullptr,                    // GAdata_ref_incrementor data_ref_incrementor
    scoreFunction,              // GAevaluate             evaluate
    randomUniformSeedInBounds,  // GAseed                 seed
    nullptr,                    // GAadapt                adapt
    nullptr,                    // GAselect_one           select_one
    nullptr,                    // GAselect_two           select_two
    nullptr,                    // GAmutate               mutate
    nullptr,                    // GAcrossover            crossover
    nullptr,                    // GAreplace              replace
    customData                  // vpointer               User data
  );

  ga_population_set_allele_min_double(pop, -DBL_MAX);

  ga_population_set_differentialevolution_parameters(
    pop,
    strategyType,
    crossoverType,
    numPerturbed,
    weightMin,
    weightMax,
    crossoverFactor
  );

  clock_gettime(CLOCK_MONOTONIC_RAW, &start);

  ga_differentialevolution(
    pop,
    generations
  );
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);

  unsigned long long elapsedTime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

  best = ga_get_entity_from_rank(pop, 0);

  const double *chromo = ((double **)best->chromosome)[0];
  printf("\n@@@\n"); // separator of run and results section
  printf("$fitness:%g$time[uS]:%llu$fitnessEvaluations:%d$generationsNeeded:%d$best:", 
          best->fitness, elapsedTime, kFitnessCalls, kGenerationsNeeded - 1);
  for (int i = 0; i < kChromoLen; i++)
  {
    printf("%g", chromo[i]);
    if (i < kChromoLen - 1)
    {
      printf(",");
    }
  }
  printf("$\n");
  ga_extinction(pop); // clean-up
  return EXIT_SUCCESS;
}
// END OF MAIN ***********************************************************************