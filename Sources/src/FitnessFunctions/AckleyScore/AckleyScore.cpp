/**
 * @file        AckleyScore.cpp
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file of Acleys function fitness score real value optimization.
 * 
 * @version     1.0
 * 
 * @date        2019-10-20 (created) \n
 *              2021-02-16 (revised)
 */

#include <FitnessFunctions/FitnessFunction.hpp>
#include <FitnessFunctions/AckleyScore.hpp>

#include <vector>
#include <cmath>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

#ifndef M_E
  #define M_E 2.71828174591064453125
#endif

AckleyScore::AckleyScore(const unsigned problemDimension,
                         const float optValue,
                         const float epsilon):
  FitnessFunction(optValue, epsilon),
  mProblemDimension(problemDimension)
{
  for (int i = 0; i < problemDimension; i++)
  {
    mUpperBounds.push_back(30.0f);
    mLowerBounds.push_back(-15.0f);
  }
}

/**
 * @brief calculates the fitness score of a given chromosome by Acleys function
 *        optimum is at x_i = 0 forall genes, domain defined for all x_i in [-15, 30]
 * 
 * @param [in, out]     chromosome       - array of size n of genes to evaluate
 * @return float                         - fitness score of given chromosome
 */
float AckleyScore::operator()(const std::vector<float>& chromosome)
{
  const float gamma             = 2.0 * M_PI;
  float       squaredInputsSum  = 0.0;
  float       cosSum            = 0.0;

  for (unsigned i = 0; i < chromosome.size(); i++)
  {
      float value       = chromosome[i];
      squaredInputsSum  += value * value;
      cosSum            += std::cos(gamma * value);
  }

  float inputAvgCoeff = 1.0f / static_cast<float>(chromosome.size());
  return (20.0f - 20.0f * std::exp(-0.2f * std::sqrt(inputAvgCoeff * squaredInputsSum)) + M_E - std::exp(inputAvgCoeff * cosSum));
}
// end of fitnessFunction **********************************************************

/**
 * @brief Way of handling the current variables that are outside of bounds
 * 
 * @param [in, out] data  - variables inside the search space to constraint.
 */
void AckleyScore::applyConstraints(std::vector<float>& data)
{
  float value = 0.0;
  for (int i = 0; i < data.size(); i++)
  {
    value = data[i];
    if (value < mLowerBounds[i])
    {
      data[i] = mLowerBounds[i];
    }
    else if (value > mUpperBounds[i])
    {
      data[i] = mUpperBounds[i];
    }
  }
}
//end of applyConstraints **********************************************************

/**
 * @brief Way for the user to get the constraints on searched variables (domain constraints)
 * 
 * @param [in] n                  - dimension of the problem (number of searched variables)
 * @param [in, out] upperBounds   - upper constraints on search space for each variable
 * @param [in, out] lowerBounds   - lower constraints on search space for each variable
 */
void AckleyScore::getConstraints(const unsigned n, 
                    std::vector<float>& upperBounds,  // out var
                    std::vector<float>& lowerBounds)  // out var
{
  for (int i = 0; i < n; i++)
  {
    upperBounds[i] = mUpperBounds[i];
    lowerBounds[i] = mLowerBounds[i];
  }
}
//end of getConstraints **********************************************************

/**
 * @brief Checks if the given chromozome is inside the constraints without modifying it.
 * 
 * @param chromozome    - the chromozome data 
 * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
 */
int AckleyScore::isInConstraints(const std::vector<float>& chromozome)
{

  int retVal = -1;
  for (int i = 0; i < chromozome.size(); i)
  {
    if (chromozome[i] > mUpperBounds[i] || chromozome[i] < mLowerBounds[i])
    {
      retVal = i;
      break;
    }
  }
  return retVal;
}
//end of isInConstraints **********************************************************
