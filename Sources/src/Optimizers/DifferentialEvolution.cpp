//project includes
#include <Optimizers/Optimizer.hpp>
#include <Optimizers/DifferentialEvolution.hpp>

#include <DataModels/DifferentialEvolutionOptions.hpp>
#include <DataModels/OptResult.hpp>

#include <FitnessFunctions/FitnessFunction.hpp>

#include <Logging/EvoLogGenerator.hpp>
#include <Logging/LogStrings.hpp>

#include "GaulLayer.hpp"

//thirdparty
extern "C"
{
  #include <gaul/gaul.h>
}
#include <fmt/format.h>

//std
#include <vector>
#include <map>
#include <ctime>

DifferentialEvolution::DifferentialEvolution(FitnessFunction&              fitness, 
                                             DifferentialEvolutionOptions& difEvoOptions,
                                             OutputFormat                  logFormat):
  Optimizer(fitness, difEvoOptions.mGeneralOptions, logFormat),
  mDifferentialEvolutionOptions(difEvoOptions)
{
  init();
}

DifferentialEvolution::~DifferentialEvolution()
{
  ga_extinction(pop); // clean-up
  delete[] customData;
}

OptResult DifferentialEvolution::run()
{
  mRunParams.mStartTimeStamp = std::time(nullptr);

  ga_differentialevolution(
    pop,
    mRunParams.mMaxGenerations
  );

  finalize();

  return OptResult
  {
    mRunParams.mBestEver,
    mRunParams.mBestEverChromosome
  };
}

void DifferentialEvolution::init() 
{
  customData = new void*[1];
  customData[0]     = this;
  //customData[1]     = new std::bind(&DifferentialEvolution::printProgress, this);

  random_seed(mRunParams.mSeed);

  pop = ga_genesis_double(
    mRunParams.mPopSize,                                    // const int              population_size
    1,                                                      // const int              num_chromo
    mRunParams.mProblemDimension,                           // const int              len_chromo
    GaulLayer::gaulGenerationCallback,                      // GAgeneration_hook      generation_hook
    nullptr,                                                // GAiteration_hook       iteration_hook
    nullptr,                                                // GAdata_destructor      data_destructor
    nullptr,                                                // GAdata_ref_incrementor data_ref_incrementor
    GaulLayer::gaulScoreCallback,                           // GAevaluate             evaluate
    GaulLayer::gaulUniformBoundedSeedCallback,              // GAseed                 seed
    nullptr,                                                // GAadapt                adapt
    nullptr,                                                // GAselect_one           select_one
    nullptr,                                                // GAselect_two           select_two
    nullptr,                                                // GAmutate               mutate
    nullptr,                                                // GAcrossover            crossover
    nullptr,                                                // GAreplace              replace
    customData                                              // vpointer               User data
  );

  ga_population_set_allele_min_double(pop, -std::numeric_limits<double>::max());

  ga_population_set_differentialevolution_parameters(
    pop,
    static_cast<ga_de_strategy_type>( mDifferentialEvolutionOptions.mStrategy),
    static_cast<ga_de_crossover_type>(mDifferentialEvolutionOptions.mCrossoverType),
    mDifferentialEvolutionOptions.mNumPerturbed,
    mDifferentialEvolutionOptions.mWeightingMin,
    mDifferentialEvolutionOptions.mWeightingMax,
    mDifferentialEvolutionOptions.mCrossoverFactor
  );

  std::map<std::string, std::string> customLogData
  {
    { LogStrings::DEString::mStrategyTypeString,     getStrategyTypeName(  mDifferentialEvolutionOptions.mStrategy)        },
    { LogStrings::DEString::mCrossoverTypeString,    getCrossoverTypeName( mDifferentialEvolutionOptions.mCrossoverType)   },
    { LogStrings::DEString::mNumPerturbedString,     fmt::format("{d}",    mDifferentialEvolutionOptions.mNumPerturbed)    },
    { LogStrings::DEString::mWeightingMinString,     fmt::format("{g}",    mDifferentialEvolutionOptions.mWeightingMin)    },
    { LogStrings::DEString::mWeightingMaxString,     fmt::format("{g}",    mDifferentialEvolutionOptions.mWeightingMax)    },
    { LogStrings::DEString::mCrossoverFactorString,  fmt::format("{g}",    mDifferentialEvolutionOptions.mCrossoverFactor) }
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

void DifferentialEvolution::finalize() 
{
  //float best = GaulLayer::bestEver;
  entity* best = ga_get_entity_from_rank(pop, 0);

  for (int i = 0; i < mRunParams.mProblemDimension; i++)
  {
    mRunParams.mBestEverChromosome[i] = static_cast<float>(reinterpret_cast<double**>(best->chromosome)[0][i]);
  }

  std::time_t elapsedSeconds  = std::time(nullptr) - mRunParams.mStartTimeStamp;
  std::string resultLog       = mLogGenerator.createResultsLogString(-1.0f * best->fitness, 
                                                                      elapsedSeconds, 
                                                                      mRunParams.mFitnessCallsPerformed, 
                                                                      mRunParams.mGenerationsPerformed,
                                                                      mRunParams.mBestEverChromosome);
  fmt::print(mRunParams.mLogFile, "{:s}", resultLog);
}

/**
 * @brief Function to map a Differential Evolution strategy type number to more readable format.
 * @param [in]   strat          - selection strategy type number of GAUL Differential Evolution
 * @return std::string          - a string describing Differential Evolution selection strategy type
 */
const std::string DifferentialEvolution::getStrategyTypeName(DifferentialEvolutionOptions::DifferentialEvolutionStrategy strat) const
{
  switch (strat)
  {
    case DifferentialEvolutionOptions::DifferentialEvolutionStrategy::Best:        return "Best";
    case DifferentialEvolutionOptions::DifferentialEvolutionStrategy::Rand:        return "Rand";
    case DifferentialEvolutionOptions::DifferentialEvolutionStrategy::RandToBest:  return "RandToBest";
    default: return "Undefined";
  }
}
// end of getStrategyTypeByName *********************************************************************************

/**
 * @brief Function to map a Differential Evolution crossover type number to more readable format.
 * @param [in]   cross          - crossover type number of GAUL Differential Evolution
 * @return std::string          - a string describing Differential Evolution crossover type
 */
const std::string DifferentialEvolution::getCrossoverTypeName(DifferentialEvolutionOptions::DifferentialEvolutionCrossover cross) const
{
  switch (cross)
  {
    case DifferentialEvolutionOptions::DifferentialEvolutionCrossover::Binomial:     return "Binomial";
    case DifferentialEvolutionOptions::DifferentialEvolutionCrossover::Exponential:  return "Exponential";
    default: return "Undefined";
  }
}
// end of getCrossoverTypeByName *********************************************************************************

float DifferentialEvolution::score(std::vector<float>& x) {}
float DifferentialEvolution::score(float* x, const size_t ptrLen) {}
void  DifferentialEvolution::printProgress() {}