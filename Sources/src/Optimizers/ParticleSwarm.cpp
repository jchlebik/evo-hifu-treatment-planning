//project includes
#include <Optimizers/Optimizer.hpp>
#include <Optimizers/ParticleSwarm.hpp>

#include <DataModels/ParticleSwarmOptions.hpp>
#include <DataModels/OptResult.hpp>

#include <FitnessFunctions/FitnessFunction.hpp>

#include <Logging/EvoLogGenerator.hpp>
#include <Logging/LogStrings.hpp>

//thirdparty
#include <libopt/opt.h> // extern C already declared inside for some reason
extern "C"
{
  #include <libopt/common.h>
  #include <libopt/pso.h>
}
#include <fmt/format.h>

//std
#include <vector>
#include <map>
#include <ctime>

ParticleSwarm::ParticleSwarm(FitnessFunction& fitness, 
                             ParticleSwarmOptions& psoOptions,
                             OutputFormat logFormat):
  Optimizer(fitness, psoOptions.mGeneralOptions, logFormat),
  mParticleSwarmOptions(psoOptions)
{
  init();
}

ParticleSwarm::~ParticleSwarm()
{
  auto* s = mSearchSpace.release();
  DestroySearchSpace(&s, _PSO_);
}

void ParticleSwarm::updateLocalAndGlobalBests(const float& currentFit, const int& particleIndex)
{
  /* It updates the local best value and position */
  if (currentFit < mSearchSpace->a[particleIndex]->best_fit)  
  {
    mSearchSpace->a[particleIndex]->best_fit = currentFit;
    for (int j = 0; j < mSearchSpace->n; j++)
    {
      mSearchSpace->a[particleIndex]->xl[j] = mSearchSpace->a[particleIndex]->x[j];
    }
  }

  /* It updates the global best value and position */
  if (currentFit < mSearchSpace->gfit) 
  { 
    mSearchSpace->gfit = currentFit;
    mSearchSpace->best = particleIndex;
    //mRunParams.mBestEver = currentFit;
    for (int j = 0; j < mSearchSpace->n; j++)
    {
      mSearchSpace->g[j] = mSearchSpace->a[particleIndex]->x[j];
      //mRunParams.mBestEverChromosome[particleIndex] = mSearchSpace->a[particleIndex]->x[j];
    }
  }
}

OptResult ParticleSwarm::run()
{
  mRunParams.mStartTimeStamp = std::time(nullptr);

  std::vector<float> chromo(mRunParams.mProblemDimension);

  while(mRunParams.mGenerationsPerformed < mRunParams.mMaxGenerations)
  {
    for (int i = 0; i < mSearchSpace->m; i++) 
    {
      mRunParams.mFitnessCallsPerformed++;
      std::copy_n(mSearchSpace->a[i]->x, mRunParams.mProblemDimension, chromo.begin());
      mSearchSpace->a[i]->fit = score(chromo); /* It executes the fitness function for agent i */
      updateLocalAndGlobalBests(mSearchSpace->a[i]->fit, i);
    }
    
    for (int i = 0; i < mSearchSpace->m; i++)
    {
      UpdateParticleVelocity(mSearchSpace.get(), i);      
      UpdateParticlePosition(mSearchSpace.get(), i);
    }
  
    mRunParams.mBestEver = static_cast<float>(mSearchSpace->gfit);
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

void ParticleSwarm::init() 
{
  SearchSpace* s = CreateSearchSpace(mRunParams.mPopSize, 
                                     mRunParams.mProblemDimension, 
                                     _PSO_);
  mSearchSpace.reset(s);

  mSearchSpace->iterations = mRunParams.mMaxGenerations;
  mSearchSpace->c1 = mParticleSwarmOptions.mC1;
  mSearchSpace->c2 = mParticleSwarmOptions.mC2;
  mSearchSpace->w  = mParticleSwarmOptions.mWeight;

  auto upperBounds = mFitnessFunction.getUpperBounds();
  auto lowerBounds = mFitnessFunction.getLowerBounds();

  for(int chromo = 0; chromo < mRunParams.mProblemDimension; chromo++)
  {
    mSearchSpace->LB[chromo] = lowerBounds[chromo % lowerBounds.size()];
    mSearchSpace->UB[chromo] = upperBounds[chromo % upperBounds.size()];
  }
  InitializeSearchSpace(mSearchSpace.get(), _PSO_);

  std::map<std::string, std::string> customLogData
  {
    { LogStrings::PSOStrings::mC1String,      fmt::format("{g}", mParticleSwarmOptions.mC1)     },
    { LogStrings::PSOStrings::mC2String,      fmt::format("{g}", mParticleSwarmOptions.mC2)     },
    { LogStrings::PSOStrings::mWeightString,  fmt::format("{g}", mParticleSwarmOptions.mWeight) }
  };

  std::string initString = mLogGenerator.createLaunchLogString(mRunParams.mPopSize, 
                                                               mRunParams.mProblemDimension, 
                                                               mRunParams.mSeed, 
                                                               std::time(nullptr),
                                                               mFitnessFunction.getUpperBounds(),
                                                               mFitnessFunction.getLowerBounds(),
                                                               customLogData);
  fmt::print(mRunParams.mLogFile, "{:s}", initString);
}

void ParticleSwarm::finalize() 
{
  for (int i = 0; i < mRunParams.mProblemDimension; i++)
  {
    mRunParams.mBestEverChromosome[i] = static_cast<float>(mSearchSpace->g[i]);
  }

  std::time_t elapsedSeconds  = std::time(nullptr) - mRunParams.mStartTimeStamp;
  std::string resultLog       = mLogGenerator.createResultsLogString(mRunParams.mBestEver, 
                                                                     elapsedSeconds, 
                                                                     mRunParams.mFitnessCallsPerformed, 
                                                                     mRunParams.mGenerationsPerformed,
                                                                     mRunParams.mBestEverChromosome);
  fmt::print(mRunParams.mLogFile, "{:s}", resultLog);
}

float ParticleSwarm::score(std::vector<float>& x) 
{
  std::vector<float> constrainedValues(x);
  mFitnessFunction.applyConstraints(constrainedValues);
  
  const float score = mFitnessFunction(constrainedValues);
  return score;
}

float ParticleSwarm::score(float* x, const size_t ptrLen) 
{
  std::vector<float> constrainedValues(x, x + ptrLen);
  mFitnessFunction.applyConstraints(constrainedValues);

  const float score = mFitnessFunction(constrainedValues);
  return score;
}

void  ParticleSwarm::printProgress()
{
  std::vector<float> popFitnesses(mRunParams.mPopSize);
  for (int i = 0; i < mRunParams.mPopSize; i++)
  {
    popFitnesses[i] = mSearchSpace->a[i]->fit;
  }

  std::string genString = mLogGenerator.createEvolutionLogString(mRunParams.mGenerationsPerformed, 
                                                                 popFitnesses, 
                                                                 mRunParams.mBestEver);
  fmt::print(mRunParams.mLogFile, "{:s}", genString);
}