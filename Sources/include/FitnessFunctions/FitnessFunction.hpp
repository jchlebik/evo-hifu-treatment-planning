/**
 * @file        FitnessFunction.hpp
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file for defining FitnessFunction class to derive by fitness implementations.
 * 
 * @version     1.0
 * 
 * @date        2019-10-20 (created) \n
 *              2021-02-16 (revised)
 */

#pragma once
#ifndef KEVOOPT_FITNESSFUNCTIONS_FITNESSFUNCTION_HPP_INCLUDE
#define KEVOOPT_FITNESSFUNCTIONS_FITNESSFUNCTION_HPP_INCLUDE

#include <vector>
#include <cmath>

/**
 * @class   FitnessFunction
 * @brief   An abstract class describing the interface of the fitness function. Functor.
 *
 */
class FitnessFunction
{
protected:
  /// the difference between two solutions
  const float         mEpsilon;
  /// known optimal value for given problem
  const float         mOptimum; 
  /// Upper bounds of the search space.
  std::vector<float>  mUpperBounds;
  /// Lower bounds of the search space.
  std::vector<float>  mLowerBounds;

public:
  FitnessFunction(const float optValue,
                  const float epsilon): 
    mOptimum(optValue),
    mEpsilon(epsilon)
  {}

  FitnessFunction(const FitnessFunction& other):
    mEpsilon(other.mEpsilon),
    mOptimum(other.mOptimum),
    mUpperBounds(other.mUpperBounds),
    mLowerBounds(other.mLowerBounds)
  {}

  inline const std::vector<float>& getUpperBounds()
  {
    return mUpperBounds;
  }

  inline const std::vector<float>& getLowerBounds()
  {
    return mLowerBounds;
  }

  /**
   * @brief calculates the fitness score of a given chromosome by the implemented function.
   * 
   * @param [in]     chromosome       - vector of genes to evaluate.
   * @return float                    - fitness score of given chromosome.
   */
  virtual float operator()(const std::vector<float>& chromosome) = 0;

  /**
   * @brief Checks the distance of current best solution from the optimum.
   * 
   * @return boolean              - was the optimum found.
   */
  virtual bool isConverged(const float& currentBest) const
  {
    return std::abs(mOptimum - currentBest) < mEpsilon;
  }

  /**
   * @brief Way of handling the current variables that are outside of bounds
   * 
   * @param [in, out] data  - variables inside the search space to constraint.
   */
  virtual void applyConstraints(std::vector<float>& data) = 0;

  // /**
  //  * @brief Way for the user to get the constraints on searched variables (domain constraints)
  //  * 
  //  * @param [in] n                  - dimension of the problem (number of searched variables)
  //  * @param [in, out] upperBounds   - upper constraints on search space for each variable
  //  * @param [in, out] lowerBounds   - lower constraints on search space for each variable
  //  */
  // virtual void getConstraints(const unsigned n, 
  //                             std::vector<float>& upperBounds,        // out var
  //                             std::vector<float>& lowerBounds) = 0;   // out var

  /**
   * @brief Checks if the given chromozome is inside the constraints without modifying it.
   * 
   * @param chromozome    - the chromozome data 
   * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
   */
  virtual int isInConstraints(const std::vector<float>& chromozome) = 0;
};

#endif //FITNESS_FUNCTION_HPP_INCLUDED