#pragma once
#ifndef HIFUEVO_INCLUDE_DATAMODELS_TABUSEARCHOPTIONS_HPP_INCLUDE
#define HIFUEVO_INCLUDE_DATAMODELS_TABUSEARCHOPTIONS_HPP_INCLUDE

#include <DataModels/GeneralOptions.hpp>

#include <cstdio>

class TabuSearchOptions
{
public:
  GeneralOptions mGeneralOptions;

  const unsigned int mSearchCount;
  const unsigned int mListLength;
  const float        mPerAllelMutationRate;
  const float        mMutationStdDev;

  TabuSearchOptions() = delete;

  TabuSearchOptions(const unsigned int& searchCount,
                    const unsigned int& listLength,
                    const float& perAllelMutationRate,
                    const float& mutationStdDev):
    mSearchCount(searchCount),
    mListLength(listLength),
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