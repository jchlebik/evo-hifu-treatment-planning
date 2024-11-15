/**
 * @file        CMAES.cpp
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for CMAES solver
 * 
 * @version     1.2
 * 
 * @date        2019-01-09 (created) \n
 *              2021-03-14 (revised)
 */

//project includes
#include <Optimizers/Optimizer.hpp>
#include <Optimizers/CovarianceMatrixAdaptationES.hpp>

#include <DataModels/CovarianceMatrixAdaptationESOptions.hpp>
#include <DataModels/OptResult.hpp>

#include <FitnessFunctions/FitnessFunction.hpp>

#include <Logging/EvoLogGenerator.hpp>

//thirdparty
#include <cmaespp/cmaes.h>
#include <fmt/format.h>

//std
#include <vector>
#include <string>
#include <algorithm>
#include <random>

#include <cstdlib>
#include <ctime>
#include <cmath>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

CovarianceMatrixAdaptationES::CovarianceMatrixAdaptationES(FitnessFunction& fitness, 
                                                           CovarianceMatrixAdaptationESOptions& cmaesOptions,
                                                           OutputFormat logFormat):
  Optimizer(fitness, cmaesOptions.mGeneralOptions, logFormat),
  mCovarianceMatrixAdaptationESOptions(cmaesOptions)
{
  init();
}

CovarianceMatrixAdaptationES::~CovarianceMatrixAdaptationES()
{
  mEvo->~CMAES();
  mCmaesParams->~Parameters();
  
  mEvo.release();
  mCmaesParams.release();
}

// void CovarianceMatrixAdaptationES::setPopSize(const unsigned newPopSize)
// {
//   mRunParams.mPopSize = newPopSize;
  
//   mCovarianceMatrixAdaptationESOptions.mStdDevs.resize(mRunParams.mPopSize);
//   mCovarianceMatrixAdaptationESOptions.mXStart.resize(mRunParams.mPopSize);
  
//   std::uniform_real_distribution<float> distribution(0.1f, 10.0f);

//   //normalizing the range of every axis to 0 to 10
//   //later conversion by using the formula:
//   //lowerBound[i] + (upperBound[i] - lowerBound[i]) * (1.0 - std::cos(M_PI * x[i]/10.0))/2.0;
//   for(int i = 0; i < mRunParams.mProblemDimension; i++)
//   {
//     mCovarianceMatrixAdaptationESOptions.mXStart[i]   = distribution(mRunParams.mRngEngine); // start at random position between bounds
//     mCovarianceMatrixAdaptationESOptions.mStdDevs[i]  = 10.0f/3.0f;
//   }

//   mCmaesParams->~Parameters();
//   mCmaesParams.reset(new CMAESpp::Parameters<float>());

//   mCmaesParams->lambda = newPopSize;
//   mCmaesParams->init(mRunParams.mProblemDimension, 
//                      mCovarianceMatrixAdaptationESOptions.mXStart.data(), 
//                      mCovarianceMatrixAdaptationESOptions.mStdDevs.data());
//   mCmaesParams->stopMaxIter = mRunParams.mMaxGenerations;

//   mEvo->~CMAES();
//   mEvo.reset(new CMAESpp::CMAES<float>());

//   mEvo->init(*mCmaesParams, mRunParams.mSeed); // freed by cmaes destructor
// }

OptResult CovarianceMatrixAdaptationES::run()
{
  float* const*       pop = nullptr;
  std::vector<float>  arFunvals(mRunParams.mPopSize);

  mRunParams.mStartTimeStamp = std::time(nullptr);

  while (mRunParams.mGenerationsPerformed < mRunParams.mMaxGenerations)
  {
    // Generate lambda new search points, sample population
    pop = mEvo->samplePopulation();
    for (int i = 0; i < mRunParams.mPopSize; i++)
    {
      mRunParams.mFitnessCallsPerformed++;
      arFunvals[i] = score(pop[i], mRunParams.mProblemDimension);
    }
    mEvo->updateDistribution(arFunvals.data());

    mRunParams.mBestEver = mEvo->get(CMAESpp::CMAES<float>::FBestEver);
    printProgress();
    mRunParams.mGenerationsPerformed++;

    if (stoppingCriteriaMet(mRunParams.mBestEver)) // if we ended set the variable
    {
      break;
    }
  }

  finalize();

  return OptResult 
  { 
    mRunParams.mBestEver,
    mRunParams.mBestEverChromosome
  };
}

void CovarianceMatrixAdaptationES::init()
{
  mCovarianceMatrixAdaptationESOptions.mStdDevs.resize(mRunParams.mProblemDimension, 0.0f);
  mCovarianceMatrixAdaptationESOptions.mXStart.resize(mRunParams.mProblemDimension, 0.0f);

  std::uniform_real_distribution<float> distribution(0.1f, 10.0f);

  //normalizing the range of every axis to 0 to 10
  //later conversion by using the formula:
  //lowerBound[i] + (upperBound[i] - lowerBound[i]) * (1.0 - std::cos(M_PI * x[i]/10.0))/2.0;
  for(int i = 0; i < mRunParams.mProblemDimension; i++)
  {
    mCovarianceMatrixAdaptationESOptions.mXStart[i]   = distribution(mRunParams.mRngEngine); // start at random position between bounds
    mCovarianceMatrixAdaptationESOptions.mStdDevs[i]  = 10.0f / 3.0f;
  }

  mCmaesParams.reset(new CMAESpp::Parameters<float>());

  if (mRunParams.mPopSize > 0)
  {
    mCmaesParams->lambda = mRunParams.mPopSize;
  }
  
  mCmaesParams->init(mRunParams.mProblemDimension, 
                     mCovarianceMatrixAdaptationESOptions.mXStart.data(), 
                     mCovarianceMatrixAdaptationESOptions.mStdDevs.data());
  mCmaesParams->stopMaxIter = mRunParams.mMaxGenerations;
  mRunParams.mPopSize       = mCmaesParams->lambda;

  mEvo.reset(new CMAESpp::CMAES<float>());
  mEvo->init(*mCmaesParams, mRunParams.mSeed); // freed by cmaes destructor

  std::string initString = mLogGenerator.createLaunchLogString(mRunParams.mPopSize, 
                                                               mRunParams.mProblemDimension, 
                                                               mRunParams.mSeed, 
                                                               std::time(nullptr),
                                                               mFitnessFunction.getUpperBounds(),
                                                               mFitnessFunction.getLowerBounds());
  fmt::print(mRunParams.mLogFile, "{:s}", initString);
}

void CovarianceMatrixAdaptationES::finalize()
{
  mEvo->getInto(CMAESpp::CMAES<float>::XBestEver, mRunParams.mBestEverChromosome.size(), mRunParams.mBestEverChromosome.data());

  const std::vector<float>& upperBounds = mFitnessFunction.getUpperBounds();
  const std::vector<float>& lowerBounds = mFitnessFunction.getLowerBounds();

  for (int i = 0; i < mRunParams.mProblemDimension; i++)
  {
    float scaleFactor                   = (1.0 - std::cos(M_PI * mRunParams.mBestEverChromosome[i] * 0.1f)) * 0.5f;
    mRunParams.mBestEverChromosome[i]   = lowerBounds[i] + (upperBounds[i] - lowerBounds[i]) * scaleFactor;
  }
  std::time_t elapsedSeconds  = std::time(nullptr) - mRunParams.mStartTimeStamp;
  std::string resultLog       = mLogGenerator.createResultsLogString(mRunParams.mBestEver, 
                                                                     elapsedSeconds, 
                                                                     mRunParams.mFitnessCallsPerformed, 
                                                                     mRunParams.mGenerationsPerformed,
                                                                     mRunParams.mBestEverChromosome);
  fmt::print(mRunParams.mLogFile, "{:s}", resultLog);
}
/**
 * @brief Calculates the overall score of a given individual.
 *
 * @param [in, out] x             - the individual to evaluate.
 * @return float                  - the score given to the evaluated individual.
 */
float CovarianceMatrixAdaptationES::score(std::vector<float> & x)
{
  std::vector<float> convertedValues(mRunParams.mProblemDimension);
  int penalize = 0;

  std::vector<float> upperBounds = mFitnessFunction.getUpperBounds();
  std::vector<float> lowerBounds = mFitnessFunction.getLowerBounds();

  //covert the normalized varibale to the original axis range
  for (int i = 0; i < mRunParams.mProblemDimension; i++)
  {
    if (x[i] < 0.0 || x[i] > 10.0)
    {
      penalize++;
    }

    float scaleFactor = (1.0 - std::cos(M_PI * x[i] * 0.1f)) * 0.5f;
    convertedValues[i] = lowerBounds[i] + (upperBounds[i] - lowerBounds[i]) * scaleFactor;
  }

  const float score = mFitnessFunction(convertedValues);
  return score + score * (static_cast<float>(penalize)/static_cast<float>(mRunParams.mProblemDimension));
}

float CovarianceMatrixAdaptationES::score(float* x, 
                                          const size_t ptrLen)
{
  std::vector<float> convertedValues(mRunParams.mProblemDimension);
  int penalize = 0;

  std::vector<float> upperBounds = mFitnessFunction.getUpperBounds();
  std::vector<float> lowerBounds = mFitnessFunction.getLowerBounds();

  //covert the normalized varibale to the original axis range
  for (int i = 0; i < ptrLen; i++)
  {
    if (x[i] < 0.0 || x[i] > 10.0)
    {
      penalize++;
    }

    float scaleFactor = (1.0 - std::cos(M_PI * x[i] * 0.1f)) * 0.5f;
    convertedValues[i] = lowerBounds[i] + (upperBounds[i] - lowerBounds[i]) * scaleFactor;
  }

  const float score = mFitnessFunction(convertedValues);
  return score + score * (static_cast<float>(penalize)/static_cast<float>(mRunParams.mProblemDimension));
}


void CovarianceMatrixAdaptationES::printProgress()
{
  size_t popSize = mCmaesParams->lambda;
  std::vector<float> popFitnesses(popSize, 0.0f);


  mEvo->getInto(CMAESpp::CMAES<float>::FVals, popFitnesses.size(), popFitnesses.data());

  std::string genString = mLogGenerator.createEvolutionLogString(mRunParams.mGenerationsPerformed, 
                                                                 popFitnesses, 
                                                                 mRunParams.mBestEver);
  
  fmt::print(mRunParams.mLogFile, "{:s}", genString);
}


