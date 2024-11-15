/**
 * @file        GeneticAlgorithm.hpp
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
#ifndef HIFUEVO_OPTIMIZERS_GENETICALGORITHM_HPP_INCLUDE
#define HIFUEVO_OPTIMIZERS_GENETICALGORITHM_HPP_INCLUDE


#include <Optimizers/Optimizer.hpp>
#include <DataModels/GeneticAlgorithmOptions.hpp>

#include <vector>
#include <memory>

//Forward declarations

struct population;
class FitnessFunction;
enum class OutputFormat : int;
struct OptResult;

class GeneticAlgorithm : public Optimizer
{
  population*   pop         = nullptr;    // Population of solutions.
  void**        customData  = nullptr;

  GeneticAlgorithmOptions& mGeneticAlgorithmOptions;

public:
  GeneticAlgorithm() = delete;

  GeneticAlgorithm(FitnessFunction&         fitness, 
                   GeneticAlgorithmOptions& gaOptions,
                   OutputFormat             logFormat);

  ~GeneticAlgorithm();

  OptResult run();

protected:
  void init() ;
  void finalize();

  float score(std::vector<float>& x);
  float score(float* x, const size_t ptrLen);
  void  printProgress();

  /**
   * @brief Function to map a Genetic Algorithm elitism type number to more readable format.
   * @param [in]   elitism         - elitism type number of GAUL Genetic Algorithm
   * @return std::string           - a string describing Genetic Algorithm elitism type
   */
  const std::string getElitismTypeName(GeneticAlgorithmOptions::GeneticAlgorithmElitismType elitism) const;

};
#endif