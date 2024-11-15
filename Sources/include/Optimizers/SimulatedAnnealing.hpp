/**
 * @file        SimulatedAnnealing.hpp
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
#ifndef HIFUEVO_OPTIMIZERS_SIMULATEDANNEALING_HPP_INCLUDE
#define HIFUEVO_OPTIMIZERS_SIMULATEDANNEALING_HPP_INCLUDE

#include <Optimizers/Optimizer.hpp>
#include <DataModels/SimulatedAnnealingOptions.hpp>

#include <vector>
#include <memory>

//Forward declarations
struct population;
class FitnessFunction;
enum class OutputFormat : int;
struct OptResult;


class SimulatedAnnealing : public Optimizer
{
  population*   pop         = nullptr;    // Population of solutions.
  void**        customData  = nullptr;

  SimulatedAnnealingOptions& mSimulatedAnnealingOptions;

public:
  SimulatedAnnealing() = delete;

  SimulatedAnnealing(FitnessFunction&           fitness, 
                     SimulatedAnnealingOptions& saOptions,
                     OutputFormat               logFormat);

  ~SimulatedAnnealing();

  OptResult run();

protected:
  void init();
  void finalize();

  float score(std::vector<float>& x);
  float score(float* x, const size_t ptrLen);
  void  printProgress();
};
#endif