//project includes
#include <Optimizers/Optimizer.hpp>
#include <Optimizers/SimulatedAnnealing.hpp>

#include <DataModels/SimulatedAnnealingOptions.hpp>
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

SimulatedAnnealing::SimulatedAnnealing(FitnessFunction& fitness, 
                                       SimulatedAnnealingOptions& saOptions,
                                       OutputFormat logFormat):
  Optimizer(fitness, saOptions.mGeneralOptions, logFormat),
  mSimulatedAnnealingOptions(saOptions)
{
  init();
}

SimulatedAnnealing::~SimulatedAnnealing()
{
  ga_extinction(pop); // clean-up
  delete[] customData;
}

OptResult SimulatedAnnealing::run()
{
  mRunParams.mStartTimeStamp = std::time(nullptr);

  entity* best = ga_get_free_entity(pop);
  ga_entity_seed(pop, best);

  ga_sa(
    pop,                        // population*   pop
    best,                       // entity*       initial
    mRunParams.mMaxGenerations // const int     max_generations
  );

  finalize();

  return OptResult
  {
    mRunParams.mBestEver,
    mRunParams.mBestEverChromosome
  };
}

void SimulatedAnnealing::init() 
{
  customData        = new void*[2];
  customData[0]     = this;
  customData[1]     = &(mSimulatedAnnealingOptions.mMutationStdDev);

  random_seed(mRunParams.mSeed);

  pop = ga_genesis_double(
    0,                                              // const int              population_size
    1,                                              // const int              num_chromo
    mRunParams.mProblemDimension,                   // const int              len_chromo
    nullptr,                                        // GAgeneration_hook      generation_hook
    GaulLayer::gaulIterationCallback_SA,            // GAiteration_hook       iteration_hook
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

  ga_population_set_allele_mutation_prob(pop, mSimulatedAnnealingOptions.mPerAllelMutationRate);
  ga_population_set_allele_min_double(pop, -std::numeric_limits<double>::max());

  ga_population_set_sa_parameters(
    pop,
    ga_sa_boltzmann_acceptance,               // boltzmann or linear acceptance or custom
    mSimulatedAnnealingOptions.mInitTemp,     // double init_temp
    0.0,                                      // double final_temp
    mSimulatedAnnealingOptions.mTempStep,     // double temp_step
    -1                                        // double temp_freq ... -1 => let SA decide
  );

  std::map<std::string, std::string> customLogData
  {
    { LogStrings::SAString::mInitialTempString,   fmt::format("{g}", mSimulatedAnnealingOptions.mInitTemp)  },
    { LogStrings::SAString::mTempStepString,      fmt::format("{g}", mSimulatedAnnealingOptions.mTempStep)  }
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

void SimulatedAnnealing::finalize() 
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

float SimulatedAnnealing::score(std::vector<float>& x) {}
float SimulatedAnnealing::score(float* x, const size_t ptrLen) {}
void  SimulatedAnnealing::printProgress() {}