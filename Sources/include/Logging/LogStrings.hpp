#pragma once
#ifndef KEVOOPT_LOGGING_LOGSTRINGS_HPP_INCLUDE
#define KEVOOPT_LOGGING_LOGSTRINGS_HPP_INCLUDE

#include <string>
#include <vector>

namespace LogStrings
{
  namespace LaunchSectionStrings
  {
    static constexpr const char mDateTimeString[]       = "dateTime";     
    static constexpr const char mPopSizeString[]        = "popSize";     // Size of the population
    static constexpr const char mDimensionString[]      = "dimension";   // Dimension of the problem
    static constexpr const char mSeedString[]           = "seed";        // Seed used
    static constexpr const char mLogModeString[]        = "logMode";     // logging mode 
    static constexpr const char mLowerBoundsString[]    = "lowerBounds";     // logging mode 
    static constexpr const char mUpperBoundsString[]    = "upperBounds";     // logging mode  

    static std::vector<std::string> mOrder{
      std::string(mDateTimeString),
      std::string(mPopSizeString),
      std::string(mDimensionString),
      std::string(mSeedString),
      std::string(mLogModeString),
      std::string(mLowerBoundsString),
      std::string(mUpperBoundsString)
    };   
  };

  namespace EvolutionSectionStrings
  {
    static constexpr const char mGenString[]     = "gen";     // "Generation"
    static constexpr const char mBestString[]    = "best";    // "Minimum fitness value so far"
    static constexpr const char mMinString[]     = "min";     // "Minimal fitness value in this generation"
    static constexpr const char mMaxString[]     = "max";     // "Maximal fitness value in this generation"
    static constexpr const char mMeanString[]    = "mean";    // "Mean fitness"
    static constexpr const char mQ1String[]      = "q1";      // "Quartile 1"
    static constexpr const char mMedianString[]  = "median";  // "Median fitness"
    static constexpr const char mQ3String[]      = "q3";      // "Quartile 3"
    static constexpr const char mIqrString[]     = "iqr";     // "Interquartile range"
    static constexpr const char mChromoString[]  = "chromo";  // "Chromosome of the best individual"
    static constexpr const char mPopDumpString[] = "popDump"; // "Dump of fitness values of the entire population"

    static constexpr const char mPopDiedSmall[] = "populationDied";

    static std::vector<std::string> mOrder{
      std::string(mGenString),
      std::string(mBestString),
      std::string(mMinString),
      std::string(mMaxString),
      std::string(mMeanString),
      std::string(mQ1String),
      std::string(mMedianString),
      std::string(mQ3String),
      std::string(mIqrString),
      std::string(mPopDumpString)
    };
  };

  namespace ResultSectionStrings
  {
    static constexpr const char mResFitnessString[]     = "fitness";     // Size of the population
    static constexpr const char mTimeElapsedString[]    = "time[s]";   // Dimension of the problem
    static constexpr const char mFitEvalstring[]        = "fitnessEvals";        // Seed used
    static constexpr const char mGensNeededString[]     = "generations";     // logging mode 
    static constexpr const char mBestString[]           = "best";     // logging mode 

    static const std::vector<std::string> mOrder{
      std::string(mResFitnessString),
      std::string(mTimeElapsedString),
      std::string(mFitEvalstring),
      std::string(mGensNeededString),
      std::string(mBestString)
    };
  };

  namespace SimpleParseStrings
  {
    static constexpr const char mDescriptorValSeparator[] = ":";
    static constexpr const char mDataSeparator[]          = "$";
    static constexpr const char mGenSeparator[]           = "\n";
    static constexpr const char mSectionSeparator[]       = "@@@";
  };

  namespace CSVStrings
  {
    static constexpr const char mCsvSeparator[]           = ";";
    static constexpr const char mArraySeparator[]         = ",";
  };

  namespace DEString
  {
    static constexpr const char mStrategyTypeString[]     = "strategyType";     
    static constexpr const char mCrossoverTypeString[]    = "crossoverType";  
    static constexpr const char mNumPerturbedString[]     = "numPerturbed";  
    static constexpr const char mWeightingMinString[]     = "weightMin";  
    static constexpr const char mWeightingMaxString[]     = "weightMax";  
    static constexpr const char mCrossoverFactorString[]  = "crossoverFactor";  
  }

  namespace GAString
  {
    static constexpr const char mElitismTypeString[]        = "elitismType";     
    static constexpr const char mTournamentSizeString[]     = "tournamentSize";  
    static constexpr const char mCrossoverRateString[]      = "crossoverRate";  
    static constexpr const char mIndividualMutRateString[]  = "individualMutation";  
    static constexpr const char mAllelMutRateString[]       = "allelMutation";
  }

  namespace SAString
  {
    static constexpr const char mInitialTempString[]    = "initialTemp";     
    static constexpr const char mTempStepString[]       = "tempStep";  
  }

  namespace TabuStrings
  {
    static constexpr const char mListLengthString[]    = "listLength";     
    static constexpr const char mSearchCountString[]   = "searchCount";  
  }

  namespace PSOStrings
  {
    static constexpr const char mC1String[]     = "C1";
    static constexpr const char mC2String[]     = "C2";
    static constexpr const char mWeightString[] = "w";
  }
};

#endif