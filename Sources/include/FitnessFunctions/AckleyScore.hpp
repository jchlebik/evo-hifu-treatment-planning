/**
 * @file        AckleyScore.hpp
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file of Ackleys function fitness score real value optimization.
 * 
 * @version     1.0
 * 
 * @date        2019-10-20 (created) \n
 *              2021-02-16 (revised)
 */
#pragma once
#ifndef KEVOOPT_FITNESSFUNCTIONS_ACKLEYSCORE_HPP_INCLUDE
#define KEVOOPT_FITNESSFUNCTIONS_ACKLEYSCORE_HPP_INCLUDE

#include <FitnessFunctions/FitnessFunction.hpp>

#include <vector>
#include <limits>

/**
 * @class   AckleyScore
 * @brief   A functor class implementing the AckleyScore fitness metric.
 *
 */
class AckleyScore : public FitnessFunction
{
  /// the dimension of the problem
  const unsigned mProblemDimension;

public:
  /**
   * @brief Constructor for the AckleyScore class.
   * 
   * @param [in] problemDimension           - the dimension of the problem.
   * @param [in] optimum                    - known optimal value for the given problem.
   * @param [in] epsilon                    - the difference between two solutions.
   * @param [in] fitnessEvaluationsLimit    - maximum number of fitness evalutions before ending the optimization.
   * 
   */
  AckleyScore(const unsigned problemDimension,
              const float optValue = 0.0f,
              const float epsilon = 0.0001f);
  
  /**
   * @brief calculates the fitness score of a given chromosome by Acleys function
   *        optimum is at x_i = 0 forall genes, domain defined for all x_i in [-15, 30]
   * 
   * @param [in, out]     chromosome       - array of size n of genes to evaluate
   * @return float                         - fitness score of given chromosome
   */
  float operator()(const std::vector<float>& chromosome);

  /**
   * @brief Way of handling the current variables that are outside of bounds
   * 
   * @param [in, out] data  - variables inside the search space to constraint.
   */
  void applyConstraints(std::vector<float>& data);

  /**
   * @brief Way for the user to get the constraints on searched variables (domain constraints)
   * 
   * @param [in] n                  - dimension of the problem (number of searched variables)
   * @param [in, out] upperBounds   - upper constraints on search space for each variable
   * @param [in, out] lowerBounds   - lower constraints on search space for each variable
   */
  void getConstraints(const unsigned n, 
                      std::vector<float>& upperBounds,   // out var
                      std::vector<float>& lowerBounds);  // out var
  
  /**
   * @brief Checks if the given chromozome is inside the constraints without modifying it.
   * 
   * @param chromozome    - the chromozome data 
   * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
   */
  int isInConstraints(const std::vector<float>& chromozome);
  
};
//end of AckleyScore
#endif