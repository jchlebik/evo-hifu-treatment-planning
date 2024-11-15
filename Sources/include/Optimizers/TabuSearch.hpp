/**
 * @file        TabuSearch.hpp
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
#ifndef HIFUEVO_OPTIMIZERS_TABUSEARCH_HPP_INCLUDE
#define HIFUEVO_OPTIMIZERS_TABUSEARCH_HPP_INCLUDE



#include <Optimizers/Optimizer.hpp>
#include <DataModels/TabuSearchOptions.hpp>

#include <vector>
#include <memory>

//Forward declarations
struct population;
class FitnessFunction;
enum class OutputFormat : int;
struct OptResult;


class TabuSearch : public Optimizer
{
  population*   pop         = nullptr;    // Population of solutions.
  void**        customData  = nullptr;

  TabuSearchOptions& mTabuSearchOptions;

public:
  TabuSearch() = delete;

  TabuSearch(FitnessFunction&   fitness, 
             TabuSearchOptions& tabuOptions,
             OutputFormat       logFormat);

  ~TabuSearch();

  OptResult run();

protected:
  void init();
  void finalize();

  float score(std::vector<float>& x);
  float score(float* x, const size_t ptrLen);
  void  printProgress();
};
#endif