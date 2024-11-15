#pragma once
#ifndef KEVOOPT_DATAMODELS_DIFFERENTIALEVOLUTIONOPTIONS_HPP_INCLUDE
#define KEVOOPT_DATAMODELS_DIFFERENTIALEVOLUTIONOPTIONS_HPP_INCLUDE

#include <DataModels/GeneralOptions.hpp>

#include <cstdio>

class DifferentialEvolutionOptions
{
public:
  enum class DifferentialEvolutionStrategy  : int { Best = 1, Rand, RandToBest, Undef = 99 };
  enum class DifferentialEvolutionCrossover : int { Binomial = 1, Exponential, Undef = 99 };

  GeneralOptions mGeneralOptions;

  DifferentialEvolutionStrategy   mStrategy;
  DifferentialEvolutionCrossover  mCrossoverType;
  const int                       mNumPerturbed;
  const float                     mWeightingMin;
  const float                     mWeightingMax;
  const float                     mCrossoverFactor;

  DifferentialEvolutionOptions() = delete;

  DifferentialEvolutionOptions(DifferentialEvolutionStrategy strategy,
                               DifferentialEvolutionCrossover crossover,
                               const int& numPerturbed,
                               const float& wMin,
                               const float& wMax,
                               const float& crossoverFactor):
    mStrategy(strategy),
    mCrossoverType(crossover),
    mNumPerturbed(numPerturbed),
    mWeightingMin(wMin),
    mWeightingMax(wMax),
    mCrossoverFactor(crossoverFactor)
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