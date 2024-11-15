#include <DataModels/DifferentialEvolutionOptions.hpp>
#include <Optimizers/DifferentialEvolution.hpp>
#include <FitnessFunctions/AckleyScore.hpp>
#include <DataModels/OptResult.hpp>

#include <iostream>

int main(int argc, const char** argv)
{
  DifferentialEvolutionOptions::DifferentialEvolutionStrategy strategy = DifferentialEvolutionOptions::DifferentialEvolutionStrategy::Best;
  DifferentialEvolutionOptions::DifferentialEvolutionCrossover crossover = DifferentialEvolutionOptions::DifferentialEvolutionCrossover::Binomial;
  const int   numPerturbed = 2;
  const float wMin = 0.1;
  const float wMax = 0.9;
  const float crossoverFactor = 0.5;

  DifferentialEvolutionOptions deOptions(strategy, crossover, numPerturbed, wMin, wMax, crossoverFactor);
  
  int seed            = 123456789;
  FILE* logFile       = stdout;
  int popSize         = 30;
  int dimension       = 40;
  int maxFitEvals     = 20000;
  int maxGenerations  = 1000;
  int maxSeconds      = 8 * 60 * 60;

  deOptions.setGeneralOptions(seed, 
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

  DifferentialEvolution de(ac, deOptions, OutputFormat::JSON);
  auto res = de.run();
  std:: cout << res.fitness << std::endl;
  return 0;
}