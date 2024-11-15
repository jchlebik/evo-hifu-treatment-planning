#include <DataModels/TabuSearchOptions.hpp>
#include <Optimizers/TabuSearch.hpp>
#include <FitnessFunctions/AckleyScore.hpp>
#include <DataModels/OptResult.hpp>

#include <iostream>

int main(int argc, const char** argv)
{
  const unsigned int searchCount = 20;
  const unsigned int listLength = 40;
  const float perAllelMutationRate = 0.2;
  const float mutationStdDev = 0.2;

  TabuSearchOptions tabuOptions(searchCount, listLength, perAllelMutationRate, mutationStdDev);
  
  int seed            = 123456789;
  FILE* logFile       = stdout;
  int popSize         = 30;
  int dimension       = 40;
  int maxFitEvals     = 20000;
  int maxGenerations  = 1000;
  int maxSeconds      = 8 * 60 * 60;

  tabuOptions.setGeneralOptions(seed, 
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

  TabuSearch tabu(ac, tabuOptions, OutputFormat::JSON);
  auto res = tabu.run();
  std:: cout << res.fitness << std::endl;
  return 0;
}