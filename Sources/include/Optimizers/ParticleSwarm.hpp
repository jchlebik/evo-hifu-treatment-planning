/**
 * @file        ParticleSwarm.hpp
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header  file describing PSO class for solving real value minimization problems
 * 
 * @version     1.2
 * 
 * @date        2020-01-09 (created) \n
 *              2021-03-04 (revised)
 */
#pragma once
#ifndef HIFUEVO_OPTIMIZERS_PARTICLESWARM_HPP_INCLUDE
#define HIFUEVO_OPTIMIZERS_PARTICLESWARM_HPP_INCLUDE


#include <Optimizers/Optimizer.hpp>
#include <DataModels/ParticleSwarmOptions.hpp>

#include <vector>
#include <memory>

class FitnessFunction;
enum class OutputFormat : int;
struct SearchSpace;
struct OptResult;

class ParticleSwarm : public Optimizer
{
  std::unique_ptr<SearchSpace> mSearchSpace;

  ParticleSwarmOptions& mParticleSwarmOptions;

public:
  ParticleSwarm() = delete;

  ParticleSwarm(FitnessFunction&      fitness, 
                ParticleSwarmOptions& psoOptions,
                OutputFormat          logFormat);

  virtual ~ParticleSwarm();

  OptResult run();

protected:

  void  init();
  void  finalize();

  /**
   * @brief Calculates the overall score of a given individual.
   *
   * @param [in, out] x             - the individual to evaluate.
   * @return float                  - the score given to the evaluated individual.
   */
  float score(std::vector<float>& x);
  float score(float* x, const size_t ptrLen);
  void  printProgress();

private:
  void updateLocalAndGlobalBests(const float& currentFit, const int& particleIndex);
};

#endif