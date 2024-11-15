#pragma once
#ifndef KEVOOPT_DATAMODELS_GENETICALGORITHMOPTIONS_HPP_INCLUDE
#define KEVOOPT_DATAMODELS_GENETICALGORITHMOPTIONS_HPP_INCLUDE

#include <DataModels/GeneralOptions.hpp>

#include <cstdio>

class GeneticAlgorithmOptions
{
public:
  enum class GeneticAlgorithmElitismType : int { RankBased = 1, KeepBestOne, NoParentsSurvive, Undef = 99 };
  
  GeneralOptions mGeneralOptions;

  GeneticAlgorithmElitismType mElitismType;
  const int                   mTournamentSelectSize;
  const float                 mCrossoverRate;
  const float                 mSelectForMutationRate;
  const float                 mPerAllelMutationRate;
  const float                 mMutationStdDev;

  GeneticAlgorithmOptions() = delete;

  GeneticAlgorithmOptions(GeneticAlgorithmElitismType elitismType,
                          const int&   tournamentSize,
                          const float& crossoverRate,
                          const float& selectForMutationRate,
                          const float& perAllelMutationRate,
                          const float& mutationStdDev):
    mElitismType(elitismType),
    mTournamentSelectSize(tournamentSize),
    mCrossoverRate(crossoverRate),
    mSelectForMutationRate(selectForMutationRate),
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