#pragma once
#ifndef HIFUEVO_OPTIMIZERS_OPTIMIZER_HPP_INCLUDE
#define HIFUEVO_OPTIMIZERS_OPTIMIZER_HPP_INCLUDE

#include <FitnessFunctions/FitnessFunction.hpp>
#include <Logging/EvoLogGenerator.hpp>
#include <DataModels/GeneralOptions.hpp>

#include <vector>

class GaulLayer;
struct OptResult;

class Optimizer
{
protected:
  FitnessFunction&  mFitnessFunction;
  EvoLogGenerator   mLogGenerator;
  GeneralOptions&   mRunParams;

  friend GaulLayer;
  
  virtual void      init() = 0;
  virtual void      finalize() = 0;

  virtual OptResult run() = 0;
  virtual float     score(std::vector<float>& x) = 0;
  virtual float     score(float* x, const size_t ptrLen) = 0;
  virtual void      printProgress() = 0;


  Optimizer(FitnessFunction& fitness, 
            GeneralOptions&  runParams,
            OutputFormat     logFormat = OutputFormat::JSON):
    mFitnessFunction(fitness),
    mRunParams(runParams),
    mLogGenerator(logFormat)
  {}

  virtual ~Optimizer() {};

  virtual const bool stoppingCriteriaMet(const float bestFitnessSoFar) const
  {
    //auto now = std::chrono::time_point_cast<std::chrono::seconds>(std::chrono::system_clock::now()).time_since_epoch().count();
    auto elapsedSeconds =  std::time(nullptr) - mRunParams.mStartTimeStamp;

    bool limit      = mRunParams.mFitnessCallsPerformed > mRunParams.mMaxFitnessEvaluations; //mFitnessFunction.reachedLimit();
    bool converged  = mFitnessFunction.isConverged(bestFitnessSoFar);//std::abs(mOptimum - currentBest) < mEpsilon; //mFitnessFunction.isConverged(bestFitnessSoFar);
    bool timesUp    = elapsedSeconds  > mRunParams.mMaxSeconds; //seconds since start
    return limit || converged || timesUp;
  }

};
#endif