/**
 * @file        FitnessFunction.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file for defining fitness functions to optimize by optimizers.
 * 
 * @version     0.2
 * 
 * @date        2019-10-20 (created) \n
 *              2020-01-11 (revised)
 */
#pragma once
#ifndef FITNESS_FUNCTION_H_INCLUDED
#define FITNESS_FUNCTION_H_INCLUDED

/**
 * @brief calculates the fitness score of a given chromosome
 * 
 * @param [in]          n                - size of the chromosome / dimension of the problem
 * @param [in, out]     chromosome       - array of size n of genes to evaluate
 * @return double                        - fitness score of given chromosome
 */
double fitnessFunction(const unsigned n, 
                       double* chromosome);

/**
 * @brief has the optimisation converged and should stop
 * 
 * @param [in] genNumber              - current generation number
 * @param [in] fitnessEvaluations     - current number of fitness executions
 * @param [in] optimum                - the optimum of the fitness function
 * @param [in] bestFitness            - best fitness found so far
 * @param [in] populationSize         - size of the current population
 * @param [in] populationFitness      - fitness values of the current population
 * @return true                       - true if the optimization has converged close enough to provided optimum and therefore should end
 * @return false                      - fals if the optimization should still continue
 */
bool isConverged(const unsigned genNumber, 
                 const unsigned fitnessEvaluations,
                 const double optimum,
                 const double bestFitness,
                 const unsigned populationSize,
                 double** populationFitness);


/**
 * @brief Checks if the given chromozome is inside the constraints.
 * 
 * @param chromoLen     - length of the chromozome to check
 * @param chromozome    - the chromozome data 
 * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
 */
int isInConstraints(const unsigned chromoLen, 
                    double* chromozome);


/**
 * @brief user defined way of how to handle the current optimized variables that are outside of bounds
 * 
 * @param [in] n          - size of the data (for problem dimension 'n' and population size 'p' its 'n * p')
 * @param [in, out] data  - current data of points in search space (variables) 
 */
void applyConstraints(const unsigned n, 
                      double** data);

/**
 * @brief way for the user to set bounds on searched variables (domain constraints)
 * 
 * @param [in] n                  - dimension of the problem (number of searched variables)
 * @param [in, out] upperBounds   - upper constraints on search space for each variable
 * @param [in, out] lowerBounds   - lower constraints on search space for each variable
 */
void getConstraints(const unsigned n, 
                    double** upperBounds,
                    double** lowerBounds);
#endif //FITNESS_FUNCTION_H_INCLUDED