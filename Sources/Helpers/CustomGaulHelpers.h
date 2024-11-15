/**
 * @file        CustomGaulHeplers.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file for custom created helpers for GAUL optimizers.
 * 
 * @version     0.2
 * 
 * @date        2019-10-20 (created) \n
 *              2020-01-11 (revised)
 */


#ifndef CUSTOM_GAUL_HELPERS_H_INCLUDED
#define CUSTOM_GAUL_HELPERS_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif
    #include <gaul.h>
#ifdef __cplusplus
}
#endif

/**
 * @brief Function to map a Genetic Algorithm scheme type number to more readable format.
 * @param [in]   scheme         - scheme type number of GAUL Genetic Algorithm
 * @return const char*          - a string describing Genetic Algorithm scheme type
 */
const char* getSchemeTypeByName(ga_scheme_type scheme);

/**
 * @brief Function to map a Genetic Algorithm elitism type number to more readable format.
 * @param [in]   scheme         - elitism type number of GAUL Genetic Algorithm
 * @return const char*          - a string describing Genetic Algorithm elitism type
 */
const char* getElitismTypeByName(ga_elitism_type elitType);

/**
 * @brief Function to map a Differential Evolution crossover type number to more readable format.
 * @param [in]   cross          - crossover type number of GAUL Differential Evolution
 * @return const char*          - a string describing Differential Evolution crossover type
 */
const char* getCrossoverTypeByName(ga_de_crossover_type cross);

/**
 * @brief Function to map a Differential Evolution strategy type number to more readable format.
 * @param [in]   strat          - selection strategy type number of GAUL Differential Evolution
 * @return const char*          - a string describing Differential Evolution selection strategy type
 */
const char* getStrategyTypeByName(ga_de_strategy_type strat);

/**
 * @brief Custom callback seeding function that works with bounds passed in pop->data.
 * 
 * @param [in]      pop         - population of the first generation
 * @param [in, out] adam        - entity to seed
 * @return success
 */
bool randomUniformSeedInBounds(population *pop, 
                               entity *adam);

/**
 * @brief Custom mutation function moving points with gaussian distribution with stddev equal to \n
 *        a third of the variable range.
 * @param [in]      pop         - gaul population structure
 * @param [in, out] father      - entity to mutate
 * @param [in, out] son         - resulting mutated entity
 * @return Nothing, data passed through parameters 
 */
void mutateDoubleMultipointCustomStddev(population *pop, 
                                        entity *father, 
                                        entity *son);


void mutateSonicationCustomStddev(population *pop, 
                                  entity *father, 
                                  entity *son);
#endif // CUSTOM_GAUL_HELPERS_H_INCLUDED