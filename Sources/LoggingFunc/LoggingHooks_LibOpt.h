/**
 * @file        LoggingHooks_LibOpt.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file of logging functions used to log progression of LibOpt library optimizers.
 * 
 * @version     0.1
 * 
 * @date        2019-02-14 (created) \n
 */
#pragma once
#ifndef LOGGING_HOOKS_LIBOPT_H_INCLUDED
#define LOGGING_HOOKS_LIBOPT_H_INCLUDED


#include "LoggingHooks.h"
#include "../LibOpt_PSO.h"


/**
 * @brief gets a statistical data from the current population 
 * 
 * @param [in, out] s             - given search space containing the population
 * @return stats*                 - structure containing all relevant statistics
 */
stats* getStats(SearchSpace* s);

/**
 * @brief creates a string representation of all fitness values in given population. 
 * 
 * @param [in, out] s             - given search space containing the population
 * @return char*                  - allocated string of information to print out
 */
char* createPopulationDump(SearchSpace* ss);

/**
 * @brief creates a string of information to log for the current generation of given population. 
 *        Used with real value population based optimizers.
 * 
 * @param [in] generation         - last evaluated generation number
 * @param [in, out] s             - given search space containing the population
 * @return char*                  - allocated string of information to print out
 */
char* createPopulationLogInfo_Double(const int generation, 
                                     SearchSpace* ss);

#endif //LOGGING_HOOKS_LIBOPT_H_INCLUDED