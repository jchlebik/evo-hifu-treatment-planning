#pragma once
#ifndef KEVOOPT_DATAMODELS_CMAEVOLUTIONSTRATEGYOPTIONS_HPP_INCLUDE
#define KEVOOPT_DATAMODELS_CMAEVOLUTIONSTRATEGYOPTIONS_HPP_INCLUDE

#include <DataModels/GeneralOptions.hpp>

#include <vector>
#include <cstdio>

class CovarianceMatrixAdaptationESOptions
{
public:
  GeneralOptions mGeneralOptions;
  
  std::vector<float> mStdDevs;
  std::vector<float> mXStart;

  CovarianceMatrixAdaptationESOptions()
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