#pragma once
#ifndef KEVOOPT_DATAMODELS_SIMULATEDANNEALNIGOPTIONS_HPP_INCLUDE
#define KEVOOPT_DATAMODELS_SIMULATEDANNEALNIGOPTIONS_HPP_INCLUDE

#include <DataModels/GeneralOptions.hpp>

#include <cstdio>

class SimulatedAnnealingOptions
{
public:
  GeneralOptions mGeneralOptions;

  const float mInitTemp;
  const float mTempStep;
  const float mPerAllelMutationRate;
  const float mMutationStdDev;

  SimulatedAnnealingOptions() = delete;

  SimulatedAnnealingOptions(const float& initTemp,
                            const float& tempStep,
                            const float& perAllelMutationRate,
                            const float& mutationStdDev):
    mInitTemp(initTemp),
    mTempStep(tempStep),
    mPerAllelMutationRate(perAllelMutationRate),
    mMutationStdDev(mutationStdDev)
  {}

  void setGeneralOptions(const int& seed, 
                         FILE* logFile,
                         const int& popSize,
                         const int& dimension,
                         const int& maxFitnessEvals,
                         const int& maxGeneralEvals,
                         const unsigned long& maxSeconds)
  {
    mGeneralOptions                         = GeneralOptions(seed, dimension, logFile);
    mGeneralOptions.mPopSize                = popSize;
    mGeneralOptions.mMaxFitnessEvaluations  = maxFitnessEvals;
    mGeneralOptions.mMaxGenerations         = maxGeneralEvals;
    mGeneralOptions.mMaxSeconds             = maxSeconds;
    mGeneralOptions.mStartTimeStamp         = std::time(nullptr);    
  }
};

#endif