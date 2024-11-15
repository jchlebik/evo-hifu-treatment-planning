#include <DataModels/SimulatedAnnealingOptions.hpp>
#include <Optimizers/SimulatedAnnealing.hpp>
#include <FitnessFunctions/AckleyScore.hpp>
#include <DataModels/OptResult.hpp>

#include <iostream>

int main(int argc, const char** argv)
{
  const float initTemp = 100;
  const float tempStep = 0.1;
  const float perAllelMutationRate = 0.2f;
  const float mutationStdDev = 0.2f;

  SimulatedAnnealingOptions saOptions(initTemp, tempStep, perAllelMutationRate, mutationStdDev);
  
  int seed            = 123456789;
  FILE* logFile       = stdout;
  int popSize         = 30;
  int dimension       = 40;
  int maxFitEvals     = 20000;
  int maxGenerations  = 1000;
  int maxSeconds      = 8 * 60 * 60;

  saOptions.setGeneralOptions(seed, 
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

  SimulatedAnnealing sa(ac, saOptions, OutputFormat::JSON);
  auto res = sa.run();
  std:: cout << res.fitness << std::endl;
  return 0;
}