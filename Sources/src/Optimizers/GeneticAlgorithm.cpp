//project includdes
#include <Optimizers/Optimizer.hpp>
#include <Optimizers/GeneticAlgorithm.hpp>

#include <DataModels/GeneticAlgorithmOptions.hpp>
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

GeneticAlgorithm::GeneticAlgorithm(FitnessFunction& fitness, 
                                   GeneticAlgorithmOptions& gaOptions,
                                   OutputFormat logFormat):
  Optimizer(fitness, gaOptions.mGeneralOptions, logFormat),
  mGeneticAlgorithmOptions(gaOptions)
{
  init();
}

GeneticAlgorithm::~GeneticAlgorithm()
{
  ga_extinction(pop); // clean-up
  delete[] customData;
}

OptResult GeneticAlgorithm::run()
{
  mRunParams.mStartTimeStamp = std::time(nullptr);

  ga_evolution(
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

void GeneticAlgorithm::init() 
{
  customData        = new void*[2];
  customData[0]     = this;
  customData[1]     = &(mGeneticAlgorithmOptions.mMutationStdDev);

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
    ga_select_one_bestof2,                                  // GAselect_one           select_one
    GaulLayer::gaulTournamentSelect,                        // GAselect_two           select_two
    GaulLayer::mutateDoubleMultipoint,                      // GAmutate               mutate
    ga_crossover_double_doublepoints,                       // GAcrossover            crossover
    nullptr,                                                // GAreplace              replace
    customData                                              // vpointer               User data
  );

  ga_population_set_allele_mutation_prob(pop, mGeneticAlgorithmOptions.mPerAllelMutationRate);
  ga_population_set_allele_min_double(pop, -std::numeric_limits<double>::max());


  ga_population_set_parameters(
    pop,                                                                  // population*           pop
    ga_scheme_type::GA_SCHEME_DARWIN,                                     // const ga_scheme_type  scheme
    static_cast<ga_elitism_type>(mGeneticAlgorithmOptions.mElitismType),  // const ga_elitism_type elitism
    mGeneticAlgorithmOptions.mCrossoverRate,                              // double                crossover
    mGeneticAlgorithmOptions.mSelectForMutationRate,                      // double                mutation
    0.0f                                                                  // double                migration
  );

  std::map<std::string, std::string> customLogData
  {
    { LogStrings::GAString::mElitismTypeString,         getElitismTypeName(mGeneticAlgorithmOptions.mElitismType)           },
    { LogStrings::GAString::mTournamentSizeString,      fmt::format("{d}", mGeneticAlgorithmOptions.mTournamentSelectSize)  },
    { LogStrings::GAString::mCrossoverRateString,       fmt::format("{g}", mGeneticAlgorithmOptions.mCrossoverRate)         },
    { LogStrings::GAString::mIndividualMutRateString,   fmt::format("{g}", mGeneticAlgorithmOptions.mMutantsInPopRate)      },
    { LogStrings::GAString::mAllelMutRateString,        fmt::format("{g}", mGeneticAlgorithmOptions.mPerAllelMutationRate)  }
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

void GeneticAlgorithm::finalize() 
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


const std::string GeneticAlgorithm::getElitismTypeName(GeneticAlgorithmOptions::GeneticAlgorithmElitismType elitism) const
{
  switch (elitism)
  {
    case GeneticAlgorithmOptions::GeneticAlgorithmElitismType::RankBased :         return "RankBasedSurvival";
    case GeneticAlgorithmOptions::GeneticAlgorithmElitismType::KeepBestOne :       return "SingleBestParentSurvives";
    case GeneticAlgorithmOptions::GeneticAlgorithmElitismType::NoParentsSurvive :  return "NoParentsSurvive";
    default: return "Undefined";
  }
}

//NOT USED------------------------------------------
float GeneticAlgorithm::score(std::vector<float>& x) 
{ 
  entity e;
  for (int i = 0; i < x.size(); i++)
  {
    static_cast<float*>(e.chromosome[0])[i] = x[i];
  }
  return GaulLayer::gaulScoreCallback(pop, &e);
}
 
float GeneticAlgorithm::score(float* x, const size_t ptrLen) 
{
  entity e;
  for (int i = 0; i < ptrLen; i++)
  {
    static_cast<float*>(e.chromosome[0])[i] = x[i];
  }
  return GaulLayer::gaulScoreCallback(pop, &e);
}

void GeneticAlgorithm::printProgress() 
{
  GaulLayer::gaulGenerationCallback(mRunParams.mGenerationsPerformed, pop);
}
//-------------------------------------------------