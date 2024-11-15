/**
 * @file        CovarianceMatrixAdaptationES.hpp
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header  file describing CMAES class for solving real value minimization problems
 * 
 * @version     1.2
 * 
 * @date        2020-01-09 (created) \n
 *              2021-03-04 (revised)
 */
#pragma once
#ifndef HIFUEVO_OPTIMIZERS_COVARIANCEMATRIXADAPTATIONES_HPP_INCLUDE
#define HIFUEVO_OPTIMIZERS_COVARIANCEMATRIXADAPTATIONES_HPP_INCLUDE

#include <Optimizers/Optimizer.hpp>
#include <DataModels/CovarianceMatrixAdaptationESOptions.hpp>
#include <DataModels/OptResult.hpp>

#include <vector>
#include <memory>


namespace CMAESpp
{
  template <typename T>
  class CMAES;

  template <typename T>
  class Parameters;
};
class FitnessFunction;
enum class OutputFormat : int;

class CovarianceMatrixAdaptationES : public Optimizer
{
  std::unique_ptr<CMAESpp::CMAES<float>>      mEvo;
  std::unique_ptr<CMAESpp::Parameters<float>> mCmaesParams;

  CovarianceMatrixAdaptationESOptions& mCovarianceMatrixAdaptationESOptions;

public:
  CovarianceMatrixAdaptationES() = delete;

  CovarianceMatrixAdaptationES(FitnessFunction& fitness, 
                               CovarianceMatrixAdaptationESOptions& cmaesOptions,
                               OutputFormat logFormat);

  virtual ~CovarianceMatrixAdaptationES();

  //virtual void      setPopSize(const unsigned newPopSize);

  virtual OptResult run();

protected:

  virtual void  init();
  virtual void  finalize();

  /**
   * @brief Calculates the overall score of a given individual.
   *
   * @param [in, out] x             - the individual to evaluate.
   * @return float                  - the score given to the evaluated individual.
   */
  virtual float score(std::vector<float>& x);
  virtual float score(float* x, const size_t ptrLen);

  virtual void  printProgress();
};


#endif