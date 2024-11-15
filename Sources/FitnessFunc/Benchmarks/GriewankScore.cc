/**
 * @file        GriewankScore.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file of Griewank function fitness score for real value optimization.
 * 
 * @version     0.2
 * 
 * @date        2019-10-20 (created) \n
 *              2020-01-11 (revised)
 */

#include "FitnessFunction.h"
#include <cmath>

/**
 * @brief calculates the fitness score of a given chromosome by Griewank function
 *        optimum is at x_i = 0 forall genes, domain defined for all x_i in [-600, 600]
 * 
 * @param [in]          n                - size of the chromosome / dimension of the problem
 * @param [in, out]     chromosome       - array of size n of genes to evaluate
 * @return double                        - fitness score of given chromosome
 */
double fitnessFunction(const unsigned n, 
                       double* chromosome)
{
  double score = 0.0;
  double squaredInputsSum = 0.0;
  double cosineProduct = 1.0;

  for (int i = 0; i < n; i++)
  {
      double value = chromosome[i];
      squaredInputsSum += std::pow(value, 2);
      cosineProduct *= std::cos(value / std::sqrt(i+1));
  }

  return (1.0/4000.0 * squaredInputsSum - cosineProduct + 1);
}
// end of fitnessFunction ***********************************************************

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
                 double** populationFitness)
{
  return std::fabs(optimum - bestFitness) <= 0.001f;
}

/**
 * @brief user defined way of how to handle the current optimized variables that are outside of bounds
 * 
 * @param [in] n          - size of the data (for problem dimension 'n' and population size 'p' its 'n * p')
 * @param [in, out] data  - current data of points in search space (variables) 
 */
void applyConstraints(const unsigned dataLen, 
                      double** data)
{
  double value = 0.0;
  for (int i = 0; i < dataLen; i++)
  {
    value = (*data)[i];
    if ( value < -60.0)
    {
      (*data)[i] = -60.0;
    }
    else if (value > 60.0)
    {
      (*data)[i] = 60.0;
    }
  }
}
//END of applyConstraints

/**
 * @brief way for the user to set bounds on searched variables (domain constraints)
 * 
 * @param [in] n                  - dimension of the problem (number of searched variables)
 * @param [in, out] upperBounds   - upper constraints on search space for each variable
 * @param [in, out] lowerBounds   - lower constraints on search space for each variable
 */
void getConstraints(const unsigned n, 
                    double** upperBounds,
                    double** lowerBounds)
  
{
  for (int i = 0; i < n; i++)
  {
    (*((double**)upperBounds))[i] = 60.0;
    (*((double**)lowerBounds))[i] = -60.0;
  }
}
//END of getConstraints

/**
 * @brief Checks if the given chromozome is inside the constraints.
 * 
 * @param chromoLen     - length of the chromozome to check
 * @param chromozome    - the chromozome data 
 * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
 */
int isInConstraints(const unsigned chromoLen, 
                    double* chromozome)
{
  double *upperBounds = (double*)malloc(chromoLen * sizeof(double));
  double *lowerBounds = (double*)malloc(chromoLen * sizeof(double));
  getConstraints(chromoLen, &upperBounds, &lowerBounds);

  int retVal = -1;
  for (int i = 0; i < chromoLen; i)
  {
    if (chromozome[i] > upperBounds[i] || chromozome[i] < lowerBounds[i])
    {
      retVal = i;
      break;
    }
  }
  free(upperBounds);
  free(lowerBounds);
  return retVal;
}
//END isInConstraints