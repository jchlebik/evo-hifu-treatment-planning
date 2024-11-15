/**
 * @file        DifferentialEvolution.hpp
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       TODO
 * 
 * @version     1.2
 * 
 * @date        2019-01-05 (created) \n
 *              2020-04-04 (revised)
 */
#pragma once
#ifndef HIFUEVO_OPTIMIZERS_DIFFERENTIALEVOLUTION_HPP_INCLUDE
#define HIFUEVO_OPTIMIZERS_DIFFERENTIALEVOLUTION_HPP_INCLUDE

#include <Optimizers/Optimizer.hpp>
#include <DataModels/DifferentialEvolutionOptions.hpp>

#include <vector>
#include <memory>

//Forward declarations

struct population;
class FitnessFunction;
enum class OutputFormat : int;
struct OptResult;


class DifferentialEvolution : public Optimizer
{
  population*   pop         = nullptr;    // Population of solutions.
  void**        customData  = nullptr;

  DifferentialEvolutionOptions& mDifferentialEvolutionOptions;

public:
  DifferentialEvolution() = delete;

  DifferentialEvolution(FitnessFunction&              fitness, 
                        DifferentialEvolutionOptions& difEvoOptions,
                        OutputFormat                  logFormat);

  ~DifferentialEvolution();

  OptResult run();

protected:
  void init() ;
  void finalize();

  float score(std::vector<float>& x);
  float score(float* x, const size_t ptrLen);
  void  printProgress();

  /**
   * @brief Function to map a Differential Evolution strategy type number to more readable format.
   * @param [in]   strat          - selection strategy type number of GAUL Differential Evolution
   * @return std::string          - a string describing Differential Evolution selection strategy type
   */
  const std::string getStrategyTypeName(DifferentialEvolutionOptions::DifferentialEvolutionStrategy strat) const;

  /**
   * @brief Function to map a Differential Evolution crossover type number to more readable format.
   * @param [in]   cross          - crossover type number of GAUL Differential Evolution
   * @return std::string          - a string describing Differential Evolution crossover type
   */
  const std::string getCrossoverTypeName(DifferentialEvolutionOptions::DifferentialEvolutionCrossover cross) const;
};
#endif