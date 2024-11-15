#pragma once
#ifndef KEVOOPT_DATAMODELS_PARTICLESWARMOPTIONS_HPP_INCLUDE
#define KEVOOPT_DATAMODELS_PARTICLESWARMOPTIONS_HPP_INCLUDE

#include <DataModels/GeneralOptions.hpp>

#include <cstdio>

class ParticleSwarmOptions
{
public:
  GeneralOptions mGeneralOptions;

  const float mC1;
  const float mC2;
  const float mWeight;

  ParticleSwarmOptions() = delete;

  ParticleSwarmOptions(const float& c1,
                       const float& c2,
                       const float& w):
    mC1(c1),
    mC2(c2),
    mWeight(w)
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