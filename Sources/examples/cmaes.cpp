#include <DataModels/CovarianceMatrixAdaptationESOptions.hpp>
#include <Optimizers/CovarianceMatrixAdaptationES.hpp>
#include <FitnessFunctions/AckleyScore.hpp>

#include <iostream>
#include <cstdio>  // FILE*

int main(int argc, const char** argv)
{
  CovarianceMatrixAdaptationESOptions cmaesOptions;

  int seed            = 123456789;
  FILE* logFile       = stdout;
  int popSize         = -1;
  int dimension       = 40;
  int maxFitEvals     = 20000;
  int maxGenerations  = 1000;
  int maxSeconds      = 8 * 60 * 60;

  cmaesOptions.setGeneralOptions(seed, 
                                 logFile,
                                 popSize,
                                 dimension,
                                 maxFitEvals,
                                 maxGenerations,
                                 maxSeconds);

  float optScore = 0.0f;

  AckleyScore ac(cmaesOptions.mGeneralOptions.mProblemDimension, 
                 optScore, 
                 std::numeric_limits<float>::epsilon());

  CovarianceMatrixAdaptationES cmaes(ac, cmaesOptions, OutputFormat::JSON);
  auto res = cmaes.run();
  
  std:: cout << res.fitness << std::endl;
  return 0;
}