//project includes
#include <Optimizers/Optimizer.hpp>
#include <Optimizers/TabuSearch.hpp>

#include <DataModels/TabuSearchOptions.hpp>
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

TabuSearch::TabuSearch(FitnessFunction&   fitness, 
                       TabuSearchOptions& tabuOptions,
                       OutputFormat       logFormat):
  Optimizer(fitness, tabuOptions.mGeneralOptions, logFormat),
  mTabuSearchOptions(tabuOptions)
{
  init();
}

TabuSearch::~TabuSearch()
{
  ga_extinction(pop); // clean-up
  delete[] customData;
}

OptResult TabuSearch::run()
{
  mRunParams.mStartTimeStamp = std::time(nullptr);

  entity* best = ga_get_free_entity(pop);
  ga_entity_seed(pop, best);

  ga_tabu(
    pop,                        // population*   pop
    best,                       // entity*       initial
    mRunParams.mMaxGenerations  // const int     max_generations
  );

  finalize();

  return OptResult
  {
    mRunParams.mBestEver,
    mRunParams.mBestEverChromosome
  };
}

void TabuSearch::init() 
{
  customData        = new void*[2];
  customData[0]     = this;
  customData[1]     = &(mTabuSearchOptions.mMutationStdDev);

  random_seed(mRunParams.mSeed);

  pop = ga_genesis_double(
    0,                                              // const int              population_size
    1,                                              // const int              num_chromo
    mRunParams.mProblemDimension,                   // const int              len_chromo
    nullptr,                                        // GAgeneration_hook      generation_hook
    GaulLayer::gaulIterationCallback_Tabu,          // GAiteration_hook       iteration_hook
    nullptr,                                        // GAdata_destructor      data_destructor
    nullptr,                                        // GAdata_ref_incrementor data_ref_incrementor
    GaulLayer::gaulScoreCallback,                   // GAevaluate             evaluate
    GaulLayer::gaulUniformBoundedSeedCallback,      // GAseed                 seed
    nullptr,                                        // GAadapt                adapt
    nullptr,                                        // GAselect_one           select_one
    nullptr,                                        // GAselect_two           select_two
    GaulLayer::mutateDoubleMultipoint,              // GAmutate               mutate
    nullptr,                                        // GAcrossover            crossover
    nullptr,                                        // GAreplace              replace
    customData                                      // vpointer               User data
  );

  GaulLayer::mPopId = ga_get_population_id(pop);

  ga_population_set_allele_mutation_prob(pop, mTabuSearchOptions.mPerAllelMutationRate);
  ga_population_set_allele_min_double(pop, -std::numeric_limits<double>::max());

  ga_population_set_tabu_parameters(
    pop,
    GaulLayer::checkTabuListRoundingFloatsCallback,
    mTabuSearchOptions.mListLength,
    mTabuSearchOptions.mSearchCount
  );

  std::map<std::string, std::string> customLogData
  {
    { LogStrings::TabuStrings::mListLengthString,   fmt::format("{g}", mTabuSearchOptions.mListLength)  },
    { LogStrings::TabuStrings::mSearchCountString,  fmt::format("{g}", mTabuSearchOptions.mSearchCount) }
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

void TabuSearch::finalize() 
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

float TabuSearch::score(std::vector<float>& x) {}
float TabuSearch::score(float* x, const size_t ptrLen) {}
void  TabuSearch::printProgress() {}