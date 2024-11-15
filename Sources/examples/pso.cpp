#include <DataModels/ParticleSwarmOptions.hpp>
#include <Optimizers/ParticleSwarm.hpp>
#include <FitnessFunctions/AckleyScore.hpp>
#include <DataModels/OptResult.hpp>

#include <iostream>

int main(int argc, const char** argv)
{
  float c1  = 2.0f;
  float c2  = 2.0f;
  float w   = 0.5f;

  ParticleSwarmOptions psoOptions(c1, c2, w);

  int seed            = 123456789;
  FILE* logFile       = stdout;
  int popSize         = 30;
  int dimension       = 40;
  int maxFitEvals     = 20000;
  int maxGenerations  = 100;
  int maxSeconds      = 8 * 60 * 60;

  psoOptions.setGeneralOptions(seed, 
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

  ParticleSwarm pso(ac, psoOptions, OutputFormat::JSON);
  auto res = pso.run();

  std:: cout << res.fitness << std::endl;
  return 0;
}