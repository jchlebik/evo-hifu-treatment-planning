/**
 * @file        LoggingHooks_Gaul.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file of logging functions used to log progression of GAUL library optimizers.
 * 
 * @version     0.3
 * 
 * @date        2019-10-20 (created) \n
 *              2020-02-14 (revised)
 */

#include "LoggingHooks_Gaul.h"

/**
 * @brief gets a statistical data from the current population 
 * 
 * @param [in, out] pop           - given population to analyse
 * @return stats*                 - structure containing all relevant statistics
 */
stats* getStats(population* pop)
{
  double max, min, mean, median, variance, stddev, kurtosis, skew, q1, q3, irq;
  sort_population(pop);
  //const entity *best = ga_get_entity_from_rank(pop, 0);
  if ( pop->size < 4 || !ga_fitness_stats(pop, &max, &min, &mean, &median, &variance, &stddev, &kurtosis, &skew))
  {
    return new stats {q1 = NAN, q3 = NAN, median = NAN, mean = NAN, min = NAN, max = NAN, irq = NAN};
  }

  bool isOddLength = pop->size & 1;
  int half = pop->size / 2;
  int quarter = pop->size / 4;
  int threeFourths = pop->size * 3.0 / 4.0;
  if (isOddLength)
  {
    median = pop->entity_iarray[half]->fitness;
    q1 = (pop->entity_iarray[quarter]->fitness + pop->entity_iarray[quarter + 1]->fitness) / 2.0;
    q3 = (pop->entity_iarray[threeFourths]->fitness + pop->entity_iarray[threeFourths + 1]->fitness) / 2.0;
  }
  else
  {
    median = (pop->entity_iarray[half]->fitness + pop->entity_iarray[half + 1]->fitness) / 2.0;
    q1 = pop->entity_iarray[quarter]->fitness;
    q3 = pop->entity_iarray[threeFourths]->fitness;
  }
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
// END OF getStats ******************************************************************
/**
 * @brief creates a string representation of all entities considered  values in given population. 
 * 
 * @param [in, out] pop           - given population to print
 * @param [in, out] best          - best entity so far
 * @return char*                  - allocated string of information to print out
 */
char* createEntitiesImprovementHistoryDump(population* pop, 
                                           entity* best)
{
  static double pastBestFitness = GA_MIN_FITNESS;
  static std::string s = "";
  const char* header = "$populationDump:";
  char buffer[32];

  if (pastBestFitness != best->fitness)
  {
    pastBestFitness = best->fitness;
    snprintf(buffer, sizeof(buffer), "%g", pastBestFitness);
    s = std::string(buffer) + ((s != "") ? "," : "") + s;
  }
  
  std::string toPrint = header + s;
  char* cstr = (char*)malloc((toPrint.size() + 1)*sizeof(char));
  if (cstr == nullptr)
  {
    return nullptr;
  }
  
  strcpy(cstr, toPrint.c_str());
  return cstr;
}
// END OF createEntitiesImprovementHistoryDump ******************************************************************
/**
 * @brief creates a string representation of all entities considered  values in given population. 
 * 
 * @param [in, out] pop           - given population to print
 * @param [in, out] best          - best entity so far
 * @return char*                  - allocated string of information to print out
 */
char* createEntitiesHistoryDump(population* pop, 
                                entity* best)
{
  static double pastBestFitness = GA_MIN_FITNESS;
  static std::string s = "$populationDump:";
  char buffer[32];
  if (s[s.size()-1] != ':')
  {
    s += ",";
  }

  if (pastBestFitness != best->fitness)
  {
    pastBestFitness = best->fitness;
    snprintf(buffer, sizeof(buffer), "%g", pastBestFitness);
    s += std::string(buffer);
  }
  else
  {
    for(int i = 0; i < pop->size; i++)
    {
      if (pop->entity_iarray[i] != nullptr 
          && pop->entity_iarray[i]->fitness != GA_MIN_FITNESS 
          && pop->entity_iarray[i]->fitness != pastBestFitness)
      {
        snprintf(buffer, sizeof(buffer), "%g", pop->entity_iarray[i]->fitness);
        s += std::string(buffer);
        break;
      }
    }
  }
  
  char* cstr = (char*)malloc(s.size()+1*sizeof(char));
  if (cstr == nullptr)
  {
    return nullptr;
  }
  strcpy(cstr, s.c_str());
  return cstr;
}

// END OF createEntitiesHistoryDump ******************************************************************
/**
 * @brief creates a string representation of all fitness values in given population. 
 * 
 * @param [in, out] pop           - given population to print
 * @return char*                  - allocated string of information to print out
 */
char* createPopulationDump(population *pop)
{
  int i = 0; 
  std::string s = "$populationDump:";
  char buffer[32];

  for(int i = 0; i < pop->size; i++)
  {
    if (pop->entity_iarray[i] != nullptr && pop->entity_iarray[i]->fitness != GA_MIN_FITNESS)
    {
      snprintf(buffer, sizeof(buffer), "%g", pop->entity_iarray[i]->fitness);
      s += std::string(buffer);
      if (i + 1 != pop->size)
      {
        s += ",";
      }
    }
  }
  char* cstr = (char*)malloc(s.size()+1*sizeof(char));
  if (cstr == nullptr)
  {
    return nullptr;
  }
  strcpy(cstr, s.c_str());
  return cstr;
}
// END OF createPopulationDump ******************************************************************
/**
 * @brief creates a string of information to log for the current generation of given population. 
 *        Used with integer population based optimizers.
 * 
 * @param [in] generation         - last evaluated generation number
 * @param [in, out] pop           - current population of algorithm
 * @return char*                  - allocated string of information to print out
 */
char* createPopulationLogInfo_Integer(const int generation, 
                                            population *pop)
{
  stats* s = getStats(pop);
  if (!s)
  {
    char* res = (char*)malloc(sizeof("-1"));
    sprintf(res, "-1");
    return res;
  }
  const entity *best = ga_get_entity_from_rank(pop, 0);

  int requiredBytes = snprintf(nullptr, 0,
                               "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
                               generation, s->min, s->max, best->fitness, s->mean, s->median, s->q1, s->q3, s->irq);
  char* buffer = new(std::nothrow) char[requiredBytes+1];
  if (!buffer)
  {
      return nullptr;
  }
  sprintf(buffer, "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
          generation, s->min, s->max, best->fitness, s->mean, s->median, s->q1, s->q3, s->irq);

  std::string textToLog = appendChromoData<int>(std::string(buffer), ((int**)best->chromosome)[0], pop->len_chromosomes, "%d");
  delete s;
  char* cstr = (char*)malloc(textToLog.size()+1*sizeof(char));
  if (cstr == nullptr)
  {
    return nullptr;
  }
  strcpy(cstr, textToLog.c_str());
  return cstr;
}
// END OF createPopulationLogInfo_Integer ******************************************************************

/**
 * @brief creates a string of information to log for the current generation of given population. 
 *        Used with real value population based optimizers.
 * 
 * @param [in] generation         - last evaluated generation number
 * @param [in, out] pop           - current population of algorithm
 * @return char*                  - allocated string of information to print out
 */
char* createPopulationLogInfo_Double(const int generation, 
                                     population *pop)
{
  stats* s = getStats(pop);
  if (!s)
  {
    char* res = (char*)malloc(sizeof("-1"));
    sprintf(res, "-1");
    return res;
  }
  const entity *best = ga_get_entity_from_rank(pop, 0);

  int requiredBytes = snprintf(nullptr, 0,
                               "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
                               generation, s->min, s->max, best->fitness, s->mean, s->median, s->q1, s->q3, s->irq);
  char* buffer = new(std::nothrow) char[requiredBytes+1];
  if (!buffer)
  {
      return nullptr;
  }
  sprintf(buffer, "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
          generation, s->min, s->max, best->fitness, s->mean, s->median, s->q1, s->q3, s->irq);

  std::string textToLog = appendChromoData<double>(
    std::string(buffer), 
    ((double**)best->chromosome)[0], 
    pop->len_chromosomes, 
    "%g"
  );
  delete s;
  char* cstr = (char*)malloc(textToLog.size()+1*sizeof(char));
  if (cstr == nullptr)
  {
    return nullptr;
  }
  strcpy(cstr, textToLog.c_str());
  return cstr;
}
// END OF createPopulationLogInfo_Double *******************************************************************

/**
 * @brief creates a string of information to with information about the current iteration of optimizer. 
 *        Used with integer value non-population based optimizers.
 * 
 * @param [in] iteration           - last evaluated iteration number
 * @param [in, out] best           - so far the best found solution of the optimizer
 * @param [in, out] pop            - population of the optimizer
 * @param [in] dimension           - dimension of the solved problem
 * @return char*                   - allocated string of information to print out
 */
char* createIterationLogInfo_Integer(const int iteration, 
                                           entity *best,
                                           population *pop, 
                                           const unsigned dimension)
{
  static double max, min, mean, median, q1, q3, irq, sum2, sum3, sum4;
  double fitness = best->fitness;
  double tmp;
  double current;
  int count = iteration + 1;

  if (iteration == 0)
  {
    max = min = mean = median = fitness;
    tmp = sum2 = sum3 = sum4 = 0.0;
  }
  else
  {
    if (max < fitness)
      max = fitness;
    if (min > best->fitness)
      min = fitness;
    mean = (mean * (iteration) + fitness) / (double)count;
    median = min + (max - min) / 2.0;

    q1 = min + (median - min) / 2.0;
    q3 = median + (max - median) / 2.0;
    irq = q3 - q1;

    tmp = fitness - mean;
    sum2 += tmp * tmp;
    sum3 += tmp * tmp * tmp;
    sum4 += tmp * tmp * tmp * tmp;
  }


  for(int i = 0; i < pop->size; i++)
  {
    if (pop->entity_iarray[i] != nullptr 
        && pop->entity_iarray[i]->fitness != GA_MIN_FITNESS 
        && pop->entity_iarray[i]->fitness != fitness)
    {
      current = pop->entity_iarray[i]->fitness;
      break;
    }
  }
  int requiredBytes = snprintf(nullptr, 0,
                               "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
                               iteration, min, max, current, mean, median, q1, q3, irq);
  char* buffer = new(std::nothrow) char[requiredBytes+1];
  if (!buffer)
  {
      return nullptr;
  }
  sprintf(buffer, "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
          iteration, min, max, current, mean, median, q1, q3, irq);

  std::string textToLog = appendChromoData<int>(
    std::string(buffer), 
    ((int**)best->chromosome)[0], 
    pop->len_chromosomes, 
    "%d"
  );
  char* cstr = (char*)malloc(textToLog.size()+1*sizeof(char));
  if (cstr == nullptr)
  {
    return nullptr;
  }
  strcpy(cstr, textToLog.c_str());
  return cstr;
}
// END OF createIterationLogInfo_Integer *******************************************************************
/**
 * @brief creates a string of information to with information about the current iteration of optimizer. 
 *        Used with real value non-population based optimizers.
 * 
 * @param [in] iteration           - last evaluated iteration number
 * @param [in, out] best           - so far the best found solution of the optimizer
 * @param [in, out] pop            - population of the optimizer
 * @param [in] dimension           - dimension of the solved problem
 * @return char*                   - allocated string of information to print out
 */
char* createIterationLogInfo_Double(const int iteration, 
                                    entity *best,
                                    population* pop, 
                                    const unsigned dimension)
{
  static double max, min, mean, median, q1, q3, irq, sum2, sum3, sum4;
  double fitness = best->fitness;
  double tmp;
  double current = std::numeric_limits<double>::infinity();
  int count = iteration + 1;

  if (iteration == 0)
  {
    max = min = mean = median = fitness;
    tmp = sum2 = sum3 = sum4 = 0.0;
  }
  else
  {
    if (max < fitness)
      max = fitness;
    if (min > best->fitness)
      min = fitness;
    mean = (mean * (iteration) + fitness) / (double)count;
    median = min + (max - min) / 2.0;

    q1 = min + (median - min) / 2.0;
    q3 = median + (max - median) / 2.0;
    irq = q3 - q1;

    tmp = fitness - mean;
    sum2 += tmp * tmp;
    sum3 += tmp * tmp * tmp;
    sum4 += tmp * tmp * tmp * tmp;
  }


  for(int i = 0; i < pop->size; i++)
  {
    if (pop->entity_iarray[i] != nullptr 
        && pop->entity_iarray[i]->fitness != GA_MIN_FITNESS 
        && pop->entity_iarray[i]->fitness != fitness)
    {
      current = pop->entity_iarray[i]->fitness;
      break;
    }
  }

  if (current = std::numeric_limits<double>::infinity())
  {
    current = fitness;
  }
  int requiredBytes = snprintf(nullptr, 0,
                               "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
                               iteration, max, min, current, mean, median, q1, q3, irq);
  char* buffer = new(std::nothrow) char[requiredBytes+1];
  if (!buffer)
  {
      return nullptr;
  }
  sprintf(buffer, "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
          iteration, max, min, current, mean, median, q1, q3, irq);

  std::string textToLog = appendChromoData<double>(
    std::string(buffer), 
    ((double**)best->chromosome)[0], 
    pop->len_chromosomes, 
    "%g"
  );
  char* cstr = (char*)malloc(textToLog.size()+1*sizeof(char));
  if (cstr == nullptr)
  {
    return nullptr;
  }
  strcpy(cstr, textToLog.c_str());
  return cstr;
}
// END OF createIterationLogInfo_Double *******************************************************************