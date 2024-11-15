/**
 * @file        CustomGaulHeplers.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for custom created helpers for GAUL optimizers.
 * 
 * @version     0.2
 * 
 * @date        2019-10-20 (created) \n
 *              2020-01-11 (revised)
 */

#include "CustomGaulHelpers.h"

//**
/**
 * @brief Function to map a Differential Evolution strategy type number to more readable format.
 * @param [in]   strat          - selection strategy type number of GAUL Differential Evolution
 * @return const char*          - a string describing Differential Evolution selection strategy type
 */
const char* getStrategyTypeByName(ga_de_strategy_type strat)
{
    switch (strat)
    {
        case 0: return "GA_DE_STRATEGY_UNKNOWN";
        case 1: return "GA_DE_STRATEGY_BEST";
        case 2: return "GA_DE_STRATEGY_RAND";
        case 3: return "GA_DE_STRATEGY_RANDTOBEST";
        default: return NULL;
    }
}
// end of getStrategyTypeByName *********************************************************************************

//**
/**
 * @brief Function to map a Differential Evolution crossover type number to more readable format.
 * @param [in]   cross          - crossover type number of GAUL Differential Evolution
 * @return const char*          - a string describing Differential Evolution crossover type
 */
const char* getCrossoverTypeByName(ga_de_crossover_type cross)
{
    switch (cross)
    {
        case 0: return "GA_DE_CROSSOVER_UNKNOWN";
        case 1: return "GA_DE_CROSSOVER_BINOMIAL";
        case 2: return "GA_DE_CROSSOVER_EXPONENTIAL";
        default: return NULL;
    }
}
// end of getCrossoverTypeByName *********************************************************************************

//**
/**
 * @brief Function to map a Genetic Algorithm scheme type number to more readable format.
 * @param [in]   scheme         - scheme type number of GAUL Genetic Algorithm
 * @return const char*          - a string describing Genetic Algorithm scheme type
 */
const char* getSchemeTypeByName(ga_scheme_type scheme)
{
    switch (scheme)
    {
        case 0: return "GA_SCHEME_DARWIN";
        case 1: return "GA_SCHEME_LAMARCK_PARENTS";
        case 2: return "GA_SCHEME_LAMARCK_CHILDREN";
        case 3: return "GA_SCHEME_LAMARCK_ALL";
        case 4: return "GA_SCHEME_BALDWIN_PARENTS";
        case 8: return "GA_SCHEME_BALDWIN_CHILDREN";
        case 12: return "GA_SCHEME_BALDWIN_ALL";
        default: return NULL;
    }
}
// end of getSchemeTypeByName ***********************************************************************************

//**
/**
 * @brief Function to map a Genetic Algorithm elitism type number to more readable format.
 * @param [in]   scheme         - elitism type number of GAUL Genetic Algorithm
 * @return const char*          - a string describing Genetic Algorithm elitism type
 */
const char* getElitismTypeByName(ga_elitism_type elitType)
{
    switch (elitType)
    {
        case 0 : return "GA_ELITISM_UNKNOWN";
        case 1 : return "GA_ELITISM_PARENTS_SURVIVE";
        case 2 : return "GA_ELITISM_ONE_PARENT_SURVIVES";
        case 3 : return "GA_ELITISM_PARENTS_DIE";
        case 4 : return "GA_ELITISM_RESCORE_PARENTS";
        case 5 : return "GA_ELITISM_BEST_SET_SURVIVE";
        case 6 : return "GA_ELITISM_PARETO_SET_SURVIVE";
        default: return "UNKNOWN TYPE";
    }
}
// end of getElitismTypeByName **********************************************************************************

bool randomUniformSeedInBounds(population *pop, 
                               entity *adam)
{
  int chromo = 0;		/* Index of chromosome to seed */
  int point;		/* Index of allele to seed */

  /* Checks. */
  if (!pop) die("Null pointer to population structure passed.");
  if (!adam) die("Null pointer to entity structure passed.");

  // we pass the constraints in gaul custom data as void*
  double *lowerBounds = ((double**)pop->data)[0];
  double *upperBounds = ((double**)pop->data)[1];

  /* Seeding. */
  for (point = 0; point < pop->len_chromosomes; point++)
  {
    ((double *)adam->chromosome[chromo])[point] = random_double_range(lowerBounds[point], upperBounds[point]);
  }

  return true;
}

/**
 * @brief Custom mutation function moving points with gaussian distribution with stddev equal to \n
 *        a third of the variable range
 * @param [in]      pop         - gaul population structure
 * @param [in, out] father      - entity to mutate
 * @param [in, out] son         - resulting mutated entity
 * @return Nothing, data passed through parameters 
 */
void mutateDoubleMultipointCustomStddev(population *pop, 
                                        entity *father, 
                                        entity *son)
{
    int i;      /* Loop variable over all chromosomes */
    int chromo; /* Index of chromosome to mutate */
    int point;  /* Index of allele to mutate */

    /* Checks */
    if (!father || !son)
        die("Null pointer to entity structure passed");

    /* Copy chromosomes of parent to offspring. */
    for (i = 0; i < pop->num_chromosomes; i++)
    {
        memcpy(son->chromosome[i], father->chromosome[i], pop->len_chromosomes * sizeof(double));
    }

    double *lowerBounds = ((double**)pop->data)[0];
    double *upperBounds = ((double**)pop->data)[1];
    /*
    * Mutate by tweaking alleles.
    */
    for (chromo = 0; chromo < pop->num_chromosomes; chromo++)
    {
        for (point = 0; point < pop->len_chromosomes; point++)
        {

            if (random_boolean_prob(pop->allele_mutation_prob))
            {
                double stddev = std::abs(upperBounds[point] - lowerBounds[point])/5.0;
                ((double *)son->chromosome[chromo])[point] += random_gaussian(0, stddev);
            }
        }
    }
    return;
}
//END of mutateDoubleMultipointCustomStddev ********************************************************

/**
 * @brief Custom mutation function moving points with gaussian distribution with stddev equal to \n
 *        a third of the variable range
 * @param [in]      pop         - gaul population structure
 * @param [in, out] father      - entity to mutate
 * @param [in, out] son         - resulting mutated entity
 * @return Nothing, data passed through parameters 
 */
void mutateSonicationCustomStddev(population *pop, 
                                  entity *father, 
                                  entity *son)
{
    int i;      /* Loop variable over all chromosomes */
    int chromo; /* Index of chromosome to mutate */
    int point;  /* Index of allele to mutate */
    int sonication;

    /* Checks */
    if (!father || !son)
        die("Null pointer to entity structure passed");

    /* Copy chromosomes of parent to offspring. */
    for (i = 0; i < pop->num_chromosomes; i++)
    {
        memcpy(son->chromosome[i], father->chromosome[i], pop->len_chromosomes * sizeof(double));
    }

    double *lowerBounds = ((double**)pop->data)[0];
    double *upperBounds = ((double**)pop->data)[1];
    /*
    * Mutate by tweaking alleles.
    */
    for (chromo = 0; chromo < pop->num_chromosomes; chromo++)
    {
        for (sonication = 0; sonication < pop->len_chromosomes; sonication+=4)
        {
            if (random_boolean_prob(pop->allele_mutation_prob))
            {
                for (point = 0; point < 4; point++)
                {
                    double stddev = std::abs(upperBounds[point] - lowerBounds[point])/5.0;
                    ((double *)son->chromosome[chromo])[sonication+point] += random_gaussian(0, stddev);
                }
            }
        }
    }
    return;
}