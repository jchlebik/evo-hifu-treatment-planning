/**
 * @file        DoubleGaul_GA.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for solving real value minimization problems with GAUL library Genetic Algorithm.
 * 
 * @version     0.4
 * 
 * @date        2019-10-16 (created) \n
 *              2020-04-04 (revised)
 */

extern "C" {
  #include <gaul.h>
}

#include <time.h>
#include <sys/mman.h>

#include "FitnessFunction.h"
#include "LoggingHooks_Gaul.h"
#include "CustomGaulHelpers.h"

/**************************************************************************************
 * GLOBALS
 * We cannot pass the run arguments through data because GAUL does not provide a way to set entity data 
 * so we have to resort to global variables
**************************************************************************************/
bool    kVerboseOutput = false;   // Is the algorithm outputing data per each generation or only launch and result summaries
int     kChromoLen = 0;           // Chromozome Lenght
double  kTopScore = INFINITY;     // If known, the optimal fitness of the problem this file is trying to solve
int     kFitnessCalls;            // Number of calls to a fitness scoring function
int     kGenerationsNeeded = 0;   // Number of generations the optimizer really needed to finish
// END OF GLOBALS DECLARATION *********************************************************

/**
 * @brief calculates the overall score of a given individual.
 * 
 * @param [in, out] pop          - variable of data type representing population, defined by GAUL
 * @param [in, out] entity       - variable of data type representing currently examined individual, defined by GAUL
 * @return bool                  - should the optimization continue
 */
bool scoreFunction(population *pop, 
                   entity *entity)
{
  double *chromo = (double*)(entity->chromosome[0]);
  applyConstraints(pop->len_chromosomes, &chromo);
  double score = fitnessFunction(pop->len_chromosomes, chromo);
  kFitnessCalls += 1;
  entity->fitness = -score;
  return true;
}
// END OF scoreFunction **************************************************************

/**
 * @brief the hook function called with the start of evaluation of every generation.
 * 
 * @param [in]      generation    - a generation number currently being evaluated
 * @param [in, out] pop           - variable of data type representing population, defined by GAUL
 * @return bool                   - should the optimization continue
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
 */
void printUsage()
{
  const int elitismTypeIds[] = {0, 1, 2, 3, 4, 5, 6};
  printf("USAGE: all arguments required\n");
  printf("./... [maxGenerations] [populationSize] [problemDimension] [crossoverRate] [mutationRate] [elitismSchemeId] [optimalFitness] [rngSeed] [isVerboseOutput]\n");
  printf("Constraints for optimization set in fitness function implementing source file.\n");
  printf("Elitism types:\n");
  for (int i = 1; i < 4; i++)
  {
    printf("\t%d - %s", elitismTypeIds[i], getElitismTypeByName((ga_elitism_type)elitismTypeIds[i]));
    if (elitismTypeIds[i] == 0)
    {
      printf(" <----- DOES NOT WORK");
    }
    printf("\n");
  }
  exit(EXIT_FAILURE);
}

/**
 * @brief the main funcion of GAUL Genetic Algorithm optimizer.
 * 
 * @param [in]      argc      - a number of arguments passed to optimizer
 * @param [in]      argv      - argumenst passed to optimizer
 * @return int                - the return value of the program
 */
int main(int argc, const char **argv)
{
  population*       pop = nullptr;       // Population of solutions.
  const entity*     best = nullptr;      // Best entity found after GA.
  int               i = 0;            // Loop over alleles.
  struct timespec   start, end;

  if (argc < 2 || argc != 10 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
  {
    printUsage();
    exit(EXIT_FAILURE);
  }

  int generations = 0, popSize = 0, seed = 0;
  float crossoverRate = 0.0, mutationRate = 0.0, migrationRate = 0.0;
  ga_scheme_type schemeType = (ga_scheme_type)0;
  ga_elitism_type elitismType;

  // cmd arguments
  generations     = strtol(argv[1], nullptr, 10);
  popSize         = strtol(argv[2], nullptr, 10);
  kChromoLen      = strtol(argv[3], nullptr, 10);
  crossoverRate   = strtof(argv[4], nullptr);
  mutationRate    = strtof(argv[5], nullptr);
  //schemeType      = (ga_scheme_type)strtol(argv[6], nullptr, 10);
  elitismType     = (ga_elitism_type)strtol(argv[6], nullptr, 10);

  kTopScore       = strtof(argv[7], nullptr);
  seed            = strtol(argv[8], nullptr, 10);
  kVerboseOutput  = (bool)strtol(argv[9], nullptr, 10);

  if (elitismType < 1 || elitismType > 3 || 
      popSize < 4 || generations < 1 || 
      kChromoLen < 1 || crossoverRate < 0 || 
      crossoverRate > 1 || mutationRate < 0 || 
      mutationRate > 1)
  {
    printUsage();
    return EXIT_FAILURE;
  }

  printf("$populationSize:%d$genesInChromosome:%d$generationsNum:%d", 
          popSize, kChromoLen, generations);
  printf("$crossoverRate:%g$mutationRate:%g$evolutionSchemeId:%d$elitismSchemeId:%d$evolutionSchemeName:%s$mutationSchemeName:%s", 
          crossoverRate, mutationRate, schemeType, elitismType, 
          getSchemeTypeByName(schemeType), getElitismTypeByName(elitismType));

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
    popSize,                                // const int              population_size
    1,                                      // const int              num_chromo
    kChromoLen,                             // const int              len_chromo
    generationHook,                         // GAgeneration_hook      generation_hook
    nullptr,                                // GAiteration_hook       iteration_hook
    nullptr,                                // GAdata_destructor      data_destructor
    nullptr,                                // GAdata_ref_incrementor data_ref_incrementor
    scoreFunction,                          // GAevaluate             evaluate
    randomUniformSeedInBounds,              // GAseed                 seed
    nullptr,                                // GAadapt                adapt
    ga_select_one_bestof2,                  // GAselect_one           select_one
    ga_select_two_bestof2,                  // GAselect_two           select_two
    mutateDoubleMultipointCustomStddev,     // GAmutate               mutate
    ga_crossover_double_doublepoints,       // GAcrossover            crossover
    nullptr,                                // GAreplace              replace
    customData                              // vpointer               User data
  );

  ga_population_set_allele_mutation_prob(pop, 1.0/5.0);
  ga_population_set_allele_min_double(pop, -DBL_MAX);

  ga_population_set_parameters(
    pop,           // population*               pop
    schemeType,    // const ga_scheme_type      scheme
    elitismType,   // const ga_elitism_type     elitism
    crossoverRate, // double                    crossover
    mutationRate,  // double                    mutation
    migrationRate  // double                    migration
  );

  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  ga_evolution(
    pop,        // population*    pop
    generations // const int      max_generations
  );
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);

  unsigned long long elapsedTime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

  printf("\n@@@\n"); // separator of run and results section
  best = ga_get_entity_from_rank(pop, 0);
  const double *chromo = ((double **)best->chromosome)[0];
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