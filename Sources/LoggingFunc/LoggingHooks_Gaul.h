/**
 * @file        LoggingHooks_Gaul.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file of logging functions used to log progression of GAUL library optimizers.
 * 
 * @version     0.2
 * 
 * @date        2019-10-20 (created) \n
 *              2020-01-11 (revised)
 */
#pragma once
#ifndef LOGGING_HOOKS_GAUL_H_INCLUDED
#define LOGGING_HOOKS_GAUL_H_INCLUDED

#include "LoggingHooks.h"

#ifdef __cplusplus
extern "C" {
#endif
    #include <gaul.h>
#ifdef __cplusplus
}
#endif


/**
 * @brief gets a statistical data from the current population 
 * 
 * @param [in, out] pop           - given population to analyse
 * @return stats*                 - structure containing all relevant statistics
 */
stats* getStats(population* pop);

/**
 * @brief creates a string representation of all entities considered  values in given population. 
 * 
 * @param [in, out] pop           - given population to print
 * @param [in, out] best          - best entity so far
 * @return string                 - string of information to print out
 */
char* createEntitiesImprovementHistoryDump(population* pop, 
                                                 entity* best);

/**
 * @brief creates a string representation of all entities considered  values in given population. 
 * 
 * @param [in, out] pop           - given population to print
 * @param [in, out] best          - best entity so far
 * @return string                 - string of information to print out
 */
char* createEntitiesHistoryDump(population* pop, 
                                      entity* best);

/**
 * @brief creates a string representation of all fitness values in given population. 
 * 
 * @param [in, out] pop           - given population to print
 * @return string                 - string of information to print out
 */
char* createPopulationDump(population *pop);

/**
 * @brief creates a string of information to log for the current generation of given population. 
 *        Used with real value population based optimizers.
 * 
 * @param [in] generation         - last evaluated generation number
 * @param [in, out] pop           - current population of algorithm
 * @return char*                  - string of information to print out
 */
char* createPopulationLogInfo_Double(const int generation, 
                                     population *pop);

/**
 * @brief creates a string of information to log for the current generation of given population. 
 *        Used with integer population based optimizers.
 * 
 * @param [in] generation         - last evaluated generation number
 * @param [in, out] pop           - current population of algorithm
 * @return char*                  - string of information to print out
 */
char* createPopulationLogInfo_Integer(const int generation, 
                                            population *pop);

/**
 * @brief creates a string of information to with information about the current iteration of optimizer. 
 *        Used with real value non-population based optimizers.
 * 
 * @param [in] iteration           - last evaluated iteration number
 * @param [in, out] best           - so far the best found solution of the optimizer
 * @param [in, out] pop            - current population of algorithm
 * @param [in] dimension           - dimension of the solved problem
 * @return char*                   - string of information to print out
 */
char* createIterationLogInfo_Double(const int iteration, 
                                          entity* best, 
                                          population * pop,
                                          const unsigned dimension);

/**
 * @brief creates a string of information to with information about the current iteration of optimizer. 
 *        Used with integer value non-population based optimizers.
 * 
 * @param [in] iteration           - last evaluated iteration number
 * @param [in, out] best           - so far the best found solution of the optimizer
 * @param [in, out] pop            - current population of algorithm
 * @param [in] dimension           - dimension of the solved problem
 * @return char*                   - string of information to print out
 */
char* createIterationLogInfo_Integer(const int iteration, 
                                           entity* best, 
                                           population * pop,
                                           const unsigned dimension);

#endif //LOGGING_HOOKS_GAUL_H_INCLUDED