/**
 * @file        LoggingHooks_LibOpt.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file of logging functions used to log progression of LibOpt library optimizers.
 * 
 * @version     0.1
 * 
 * @date        2019-02-14 (created) \n
 */

//#include "LoggingHooks_LibOpt.h"

#include "LoggingHooks_LibOpt.h"

/**
 * @brief gets a statistical data from the current population 
 * 
 * @param [in, out] s             - given search space containing the population
 * @return stats*                 - structure containing all relevant statistics
 */
stats* getStats(SearchSpace* s)
{

  std::vector<Agent*> agentsVector(s->a, s->a+s->m);
  std::sort(agentsVector.begin(), agentsVector.end(), [](Agent* a1, Agent* a2){
    return a1->fit < a2->fit;
  });
  double max, min, mean = 0.0, median, /*variance, stddev, kurtosis, skew,*/ q1, q3, irq;

  const Agent *best = agentsVector[0];
  if (s->m < 4)
  {
    return new stats {q1 = NAN, q3 = NAN, median = NAN, mean = NAN, min = NAN, max = NAN, irq = NAN};
  }

  bool isOddLength = s->m & 1;
  int half = s->m / 2;
  int quarter = s->m / 4;
  int threeFourths = s->m * 3.0 / 4.0;

  if (isOddLength)
  {
    median = agentsVector[half]->fit;
    q1 = (agentsVector[quarter]->fit + agentsVector[quarter + 1]->fit) / 2.0;
    q3 = (agentsVector[threeFourths]->fit + agentsVector[threeFourths + 1]->fit) / 2.0;
  }
  else
  {
    median = (agentsVector[half]->fit + agentsVector[half + 1]->fit) / 2.0;
    q1 = agentsVector[quarter]->fit;
    q3 = agentsVector[threeFourths]->fit;
  }
  std::for_each(agentsVector.begin(), agentsVector.end(), [&] (Agent *a) {
    mean += a->fit;
  });
  mean = mean / (float)agentsVector.size();
  min = best->fit;
  max = agentsVector[s->m-1]->fit;
  irq = q3 - q1;

  stats* result = new stats();
  result->q1 = q1;
  result->q3 = q3;
  result->median = median;
  result->mean = mean;
  result->min = min;
  result->max = max;
  result->irq = irq;
  return result;
}


/**
 * @brief creates a string representation of all fitness values in given population. 
 * 
 * @param [in, out] s             - given search space containing the population
 * @return char*                  - allocated string of information to print out
 */
char* createPopulationDump(SearchSpace* ss)
{
  std::string dump = "$populationDump:";
  char buffer[32];

  for(int i = 0; i < ss->m; i++)
  {
    snprintf(buffer, sizeof(buffer), "%g", ss->a[i]->fit);
    dump += std::string(buffer);
    if (i + 1 != ss->m)
    {
      dump += ",";
    }
  }
  char* cstr = (char*)malloc(dump.size()+1*sizeof(char));
  if (cstr == nullptr)
  {
    return nullptr;
  }
  strcpy(cstr, dump.c_str());
  return cstr;
}

/**
 * @brief creates a string of information to log for the current generation of given population. 
 *        Used with real value population based optimizers.
 * 
 * @param [in] generation         - last evaluated generation number
 * @param [in, out] pop           - current population of algorithm
 * @return char*                  - allocated string of information to print out
 */
char* createPopulationLogInfo_Double(const int generation, 
                                     SearchSpace* ss)
{
  stats* st = getStats(ss);
  if (!st)
  {
    //char* res = (char*)malloc(sizeof("-1"));
    //sprintf(res, "-1");
    return nullptr;
  }
  const Agent *best = ss->a[ss->best];


  int requiredBytes = snprintf(nullptr, 0,
                               "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
                               generation, st->min, st->max, best->fit, st->mean, st->median, st->q1, st->q3, st->irq);
  char* buffer = new(std::nothrow) char[requiredBytes+1];
  if (!buffer)
  {
      return nullptr;
  }
  sprintf(buffer, "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
          generation, st->min, st->max, best->fit, st->mean, st->median, st->q1, st->q3, st->irq);

  std::string textToLog = appendChromoData<double>(
    std::string(buffer), 
    best->x, 
    ss->n, 
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