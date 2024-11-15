#include <DataModels/GeneticAlgorithmOptions.hpp>
#include <Optimizers/GeneticAlgorithm.hpp>
#include <FitnessFunctions/AckleyScore.hpp>
#include <DataModels/OptResult.hpp>

#include <iostream>

int main(int argc, const char** argv)
{
  int tournamentSize          = 2;
  float crossoverRate         = 0.8;
  float selectForMutationRate = 0.1;
  float perAllelMutationRate  = 0.3;
  float stdDev                = 0.3;

  GeneticAlgorithmOptions gaOptions(GeneticAlgorithmOptions::GeneticAlgorithmElitismType::RankBased, 
                                    tournamentSize, 
                                    crossoverRate, 
                                    selectForMutationRate, 
                                    perAllelMutationRate, 
                                    stdDev);
  int seed            = 123456789;
  FILE* logFile       = stdout;
  int popSize         = 30;
  int dimension       = 40;
  int maxFitEvals     = 20000;
  int maxGenerations  = 100;
  int maxSeconds      = 8 * 60 * 60;

  gaOptions.setGeneralOptions(seed, 
                              logFile,
                              popSize,
                              dimension,
                              maxFitEvals,
                              maxGenerations,
                              maxSeconds);
  float optScore = 0.0f;

  AckleyScore ac(dimension, 
                 optScore, 
                 std::numeric_limits<float>::epsilon());

  GeneticAlgorithm ga(ac, gaOptions, OutputFormat::JSON);
  auto res = ga.run();
  std:: cout << res.fitness << std::endl;
  return 0;
}