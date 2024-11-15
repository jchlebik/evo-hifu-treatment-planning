#pragma once
#ifndef KEVOOPT_FITNESSFUNCTIONS_HIFUSCORE_HPP_INCLUDE
#define KEVOOPT_FITNESSFUNCTIONS_HIFUSCORE_HPP_INCLUDE

#include <FitnessFunctions/FitnessFunction.hpp>
#include <FitnessFunctions/HIFUScore/DataModels/Medium.hpp> // custom dtor. unique_ptr does not like forward declares
#include <FitnessFunctions/HIFUScore/KWaveDiffusionSolver/kWaveDiffusionSolver.hpp>

#include <vector>
#include <string>
#include <limits>
#include <memory>

class DiscreteMap;
class kWaveDiffusionSolver2D;
class CircleInfo;

/**
 * @class   HIFUScore
 * @brief   A class representing the HIFU planning problem. Inherits and uses the FitnessFunction interface. Most data \n
 *          are read from a file. \n 
 *          Constraints file : \n
 *            First line in file defines   the type of constraint: 'B' for box and 'C' for circle. Second line are the \n
 *            lower bounds separated by comma. Third line are the upper bounds separated by comma. For the 'C' type, next \n 
 *            line defines the constraining circle in format of 'centerX, centerY, radius'. Bounds are periodically \n 
 *            repeated in case where the problem dimension is higher than given number of bounds. \n 
 *            Example inside 'Data/Targets/Blob/constraints.dat' \n 
 *          Medium files : \n
 *            Medium is defined by 4 files inside the Data/AustinWoman folder. \n 
 *          Target and Prohibition files: \n 
 *            Target and penalization maps are created using 2 files, one for each map. Each map is created by overlapping \n
 *            multiple circles described in format 'centerX, centerY, radius, value' inside the file. The circles are drawn \n
 *            over the medium map, where the drawn circle is filled with the given value and rest of the area is zeroed. \n
 *            For the target map, all circles are then combined and overlaps solved by function 'a > b ? a + b : b'. \n
 *            For the penalization map, all circles are combined with overlaps solved by 'a+b'. \n
 *            Examples inside 'Data/Targets/Blob/targetMapCircles.dat' and 'Data/Targets/Blob/prohibitedMapCircles.dat'
 *
 */
class HIFUScore : public FitnessFunction
{
public:
  enum  class ProblemType { Blob, Flower };

private:
  enum class ConstraintType { Box, Circle };
  /// Enumeration type for different supported hifu problem types.
  const ProblemType mProblemType;
  /// Number of sonications used to solve the problem.
  const int mNumberOfSonications;
  /// Smart pointer to original Medium class before treatment to save reinitializations.
  std::unique_ptr<Medium>                 mMediumBeforeTreatment;
  /// Smart pointer to Discrete Map class representing the target map.
  std::unique_ptr<DiscreteMap>            mTargetMap;
  /// Smart pointer to Discrete Map class representing the penalization map.
  std::unique_ptr<DiscreteMap>            mPenalizeMap;
  /// Smart pointer to kWaweDiffusionSolver2D class responsible for running the simulation.
  std::unique_ptr<kWaveDiffusionSolver2D> mKWaveDiffSolver;

  /// Enumeration of different constraints option. Box or Circle.
  ConstraintType mConstraintType;
  /// Vector of constraining circles for non-box constraints. TODO.
  std::vector<CircleInfo> mConstraintCircles;

public:
  /**
   * @brief Constructor for the HIFUScore class. Initializes and saves the original state of the medium, the target maps, \n 
   *        the diffusion solver and the constraints. As long as the reference to this instance is saved, these properties \n
   *        can and will be reused in all future evaluations. All data are read from files inside the Data folder for the \n
   *        given problem type.
   * 
   * @param [in] problemType                - the type of the problem to solve. Currently supported are Blob or Flower.
   * @param [in] numOfSonications           - number of sonications used to solve the problem.
   * @param [in] optimum                    - known optimal value for the given problem.
   * @param [in] epsilon                    - the difference between two solutions.
   * 
   */
  HIFUScore(ProblemType problemType, 
            const unsigned numOfSonications,
            const float optValue = 0.0f,
            const float epsilon = 0.0f);
  /**
   * @brief Calculates the fitness score of the given chromosome. Implemented as a functor.
   * 
   * @param [in]     chromosome       - vector of genes to evaluate
   * @return float                    - fitness score of given chromosome
   */
  float operator()(const std::vector<float>& chromosome);
  /**
   * @brief Way of handling the current variables that are outside of bounds
   * 
   * @param [in, out] data  - variables inside the search space to constraint.
   */
  void applyConstraints(std::vector<float>& data);

  /**
   * @brief Checks if the given chromozome is inside the constraints without modifying it.
   * 
   * @param chromozome    - the chromozome data 
   * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
   */
  int isInConstraints(const std::vector<float>& chromozome);

  /**
   * @brief Way for the user to get the constraints on searched variables (domain constraints)
   * 
   * @param [in] n                  - dimension of the problem (number of searched variables)
   * @param [in, out] upperBounds   - upper constraints on search space for each variable
   * @param [in, out] lowerBounds   - lower constraints on search space for each variable
   */
  void getConstraints(const unsigned n, 
                      std::vector<float>& upperBounds,    // out var
                      std::vector<float>& lowerBounds);   // out var

private:
  /**
   * @brief Way of handling the current variables that are outside of bounds. \n
   *        Specialized implementation for the Box constraitns type.
   * 
   * @param [in, out] data  - variables inside the search space to constraint.
   */
  void applyConstraintsForBox(std::vector<float>& data);

  /**
   * @brief Way of handling the current variables that are outside of bounds. \n
   *        Specialized implementation for the Circle constraitns type.
   * 
   * @param [in, out] data  - variables inside the search space to constraint.
   */
  void applyConstraintsForCircle(std::vector<float>& data);

  /**
   * @brief Checks if the given chromozome is inside the constraints without modifying it. \n 
   *        Specialized implementation for the Box constraitns type.
   * 
   * @param chromozome    - the chromozome data 
   * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
   */
  int isInConstraintsForBox(const std::vector<float>& chromozome);

  /**
   * @brief Checks if the given chromozome is inside the constraints without modifying it. \n 
   *        Specialized implementation for the Circle constraitns type.
   * 
   * @param chromozome    - the chromozome data 
   * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
   */
  int isInConstraintsForCircle(const std::vector<float>& chromozome);

  /**
   * @brief Setups the search space constraints data. All data read from a file. First line in file defines \n 
   *        the type of constraint: 'B' for box and 'C' for circle. Second line are the lower bounds separated \n 
   *        by comma. Third line are the upper bounds separated by comma. For the 'C' type, next line defines the \n
   *        constraining circle in format of 'centerX, centerY, radius'. Bounds are periodically repeated in case \n
   *        where the problem dimension is higher than given number of bounds. Example inside Data/Blob/constraints.dat
   * 
   * @param [in] constraitnsFilePath      - a file defining the constraints. 
   * 
   */
  void setupConstraints(const std::string& constraintsFilePath);

};

#endif