#pragma once
#ifndef KEVOOPT_DATAMODELS_GENERALOPTIONS_HPP_INCLUDE
#define KEVOOPT_DATAMODELS_GENERALOPTIONS_HPP_INCLUDE

#include <vector>
#include <random>

#include <ctime>   // time_t
#include <cstdio>  // FILE*

class GeneralOptions
{
public:
  //rng
  int                 mSeed;
  std::mt19937        mRngEngine;

  //general
  int                 mProblemDimension;
  int                 mPopSize;
  float               mBestEver;
  std::vector<float>  mBestEverChromosome;
  std::time_t         mStartTimeStamp;
  FILE*               mLogFile;

  //counters
  int                 mFitnessCallsPerformed;
  int                 mGenerationsPerformed;

  //stopping
  unsigned long       mMaxSeconds;
  int                 mMaxGenerations;
  int                 mMaxFitnessEvaluations;

  //distributed computing parameters
  int                 mMigrationInterval;
  int                 mNumberOfImmigrants;
  int                 mIslandId;
  int                 mIslandCount;
  int                 mTargetIsland;
  int                 mSourceIsland;
  int                 mPopSizePerIsland;

  GeneralOptions() {}

  GeneralOptions(int seed, int problemDimension, FILE* logFile = stdout) :
    mSeed(seed), 
    mRngEngine(mSeed),
    mLogFile(logFile),
    mFitnessCallsPerformed(0), 
    mGenerationsPerformed(0),
    mMaxSeconds(std::numeric_limits<unsigned long>::max()), 
    mMaxGenerations(std::numeric_limits<int>::max()), 
    mMaxFitnessEvaluations(std::numeric_limits<int>::max()),
    mStartTimeStamp(0), 
    mProblemDimension(problemDimension), 
    mPopSize(-1),
    mBestEver(std::numeric_limits<float>::max()),
    mBestEverChromosome(mProblemDimension),
    mMigrationInterval(std::numeric_limits<int>::max()),
    mNumberOfImmigrants(-1), 
    mIslandId(-1), 
    mIslandCount(-1), 
    mTargetIsland(-1), 
    mSourceIsland(-1), 
    mPopSizePerIsland(-1)
  {}
};

#endif