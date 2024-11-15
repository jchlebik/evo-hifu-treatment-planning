/**
 * @file        LoggingHooks_CMAESpp.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file of logging functions used to log progression of CMA-ESpp library.
 * 
 * @version     0.3
 * 
 * @date        2019-10-20 (created) \n
 *              2020-02-16 (revised)
 */

#include "LoggingHooks_CMAESpp.h"

/**
 * @brief creates a string of information to log for the current generation of given population. 
 *        Used with real value population based optimizers.
 * 
 * @param [in] evo           - cmaes class containing the current progress of the optimization
 * @param [in] pop           - current population
 * @return char*             - allocated string of information to print out
 */

char* createPopulationLogInfo_Double(CMAES<double> &evo,
                                    double const* pop)
{
  unsigned int generation = evo.get(CMAES<double>::Generation);
  unsigned int N = evo.get(CMAES<double>::PopSize);
  stats* st = getStats<double>(pop, N);
  if (!st)
  {
    //char* res = (char*)malloc(sizeof("-1"));
    //sprintf(res, "-1");
    return nullptr;
  }
  const double bestEver = evo.get(CMAES<double>::FBestEver);
  const double bestCurrent = evo.get(CMAES<double>::Fitness);

  const double* bestChromo = evo.getPtr(CMAES<double>::XBest);
  int requiredBytes = snprintf(nullptr, 0,
                               "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
                               generation, st->min, st->max, bestCurrent, st->mean, st->median, st->q1, st->q3, st->irq);
  char* buffer = new(std::nothrow) char[requiredBytes+1];
  if (!buffer)
  {
      return nullptr;
  }
  sprintf(buffer, "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
          generation, st->min, st->max, bestCurrent, st->mean, st->median, st->q1, st->q3, st->irq);

  std::string textToLog = appendChromoData<double>(
    std::string(buffer), 
    (double*)bestChromo, 
    (int)evo.get(CMAES<double>::Dimension), 
    "%g"
  );
  delete st;
  char* cstr = (char*)malloc((textToLog.size()+1)*sizeof(char));
  if (cstr == nullptr)
  {
    delete[] buffer;
    return nullptr;
  }
  delete[] buffer;
  strcpy(cstr, textToLog.c_str());
  return cstr;
}