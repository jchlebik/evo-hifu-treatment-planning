/**
 * @file        Logging.hpp
 *
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 *
 * @brief       The header file of logging functions used to log progression of optimizer.
 *
 * @version     1.0
 *
 * @date        2019-10-20 (created) \n
 *              2020-03-03 (revised)
 */
#pragma once
#ifndef KEVOOPT_LOGGING_EVOLOGGENERATOR_HPP_INCLUDE
#define KEVOOPT_LOGGING_EVOLOGGENERATOR_HPP_INCLUDE

#include <Logging/LogStrings.hpp>

#include <fmt/format.h>

#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <deque>
#include <iterator>
#include <ctime>    //std time
#include <utility>  //pair


using KeyValuePair = std::pair<std::string, std::string>;

enum class OutputFormat : int { CSV, JSON, SimpleParse };

class CustomLogData
{
  std::vector<KeyValuePair> mCustomData;

public:
  void add(const std::string& key, const std::string& val);

  void add(const std::pair<std::string, std::string>& stlPair);

  KeyValuePair& operator[](const int& index);

  auto begin()       ;
  auto end()         ;
  auto begin()  const;
  auto end()    const;
  auto cbegin() const;
  auto cend()   const;
  auto size()   const;
  auto empty()  const;
  
  CustomLogData(const std::map<std::string, std::string>& stlMap);

  CustomLogData(const std::vector<KeyValuePair>& stlVector);

  CustomLogData(const CustomLogData& other);

  CustomLogData(); 
};

class EvoLogGenerator
{
  bool          mCreateCSVGenHeader = true;
  OutputFormat  mOutputFormat;

  std::deque<float> mHistory;

public:
  EvoLogGenerator(OutputFormat outputFormat = OutputFormat::JSON);

  std::string createLaunchLogString(const int& popSize,
                                    const int& dimension,
                                    const long long& seed,
                                    const std::time_t dateTime,
                                    const std::vector<float>& upperBounds,
                                    const std::vector<float>& lowerBounds,                                    
                                    const CustomLogData& customData = std::map<std::string,std::string>());

  std::string createEvolutionLogString(const int&    currentGeneration,
                                       const std::vector<float>& populationFitnesses,
                                       const float&  bestSoFar,
                                       const CustomLogData& customData = std::map<std::string,std::string>());

  std::string createEvolutionLogString(const int&    currentGeneration,
                                       const float&  currentFitness,
                                       const float&  bestSoFar,
                                       const CustomLogData& customData = std::map<std::string,std::string>());

  std::string createResultsLogString(const float& bestFitness,
                                     const std::time_t& runTime,
                                     const int& fitnessEvaluations,
                                     const int& generations,
                                     const std::vector<float>& bestChromo,
                                     const CustomLogData& customData = std::map<std::string,std::string>());

private:

  std::string createCSVString(const std::vector<KeyValuePair>& pairsToPrint, const bool printHeader = false) const;

  std::string createJSONString(const std::vector<KeyValuePair>& pairsToPrint) const;

  std::string createSimpleParseString(const std::vector<KeyValuePair>& pairsToPrint) const;


  const std::map<std::string, std::string> getLaunchMap(const int& popSize,
                                                        const int& dimension,
                                                        const long long& seed,
                                                        const std::time_t dateTime,
                                                        const std::vector<float>& upperBounds,
                                                        const std::vector<float>& lowerBounds) const;

  const std::map<std::string, std::string> getGenerationMap(const int&    currentGeneration,
                                                            const std::vector<float>& populationFitnesses,
                                                            const float&  currentFitness,
                                                            const float&  bestSoFar,
                                                            const bool rolling);


  const std::map<std::string, std::string> getResultsMap(const float& bestFitness,
                                                         const std::time_t& runTime,
                                                         const int& fitnessEvaluations,
                                                         const int& generations,
                                                         const std::vector<float>& bestResult) const;


  const std::vector<KeyValuePair> orderRequiredDataAndAddCustom(const std::map<std::string, std::string> requiredData,
                                                                const std::vector<std::string>& orderingOfRequiredData,
                                                                const CustomLogData& customData) const;


  template <typename RandomIt>
  std::map<std::string, std::string> getStatisticsImpl(const RandomIt sortedBegin, 
                                                       const RandomIt sortedEnd,
                                                       std::forward_iterator_tag) const = delete;

  template <typename RandomIt>
  std::map<std::string, std::string> getStatisticsImpl(const RandomIt sortedBegin, 
                                                       const RandomIt sortedEnd,
                                                       std::bidirectional_iterator_tag) const = delete;

  /**
   * @brief gets a statistical data from the current population
   *
   * @param [in] pop            - vector of fitness values
   * @return map<string, float> - structure containing all relevant statistics
   */
  template <typename RandomIt>
  std::map<std::string, std::string> getStatisticsImpl(const RandomIt sortedBegin, 
                                                       const RandomIt sortedEnd,
                                                       std::random_access_iterator_tag) const
  {
    const int N = std::distance(sortedBegin, sortedEnd);
    float max, min, mean, median, q1, q3, iqr;

    if (N == 1)
    {
      max = min = mean = median = q1 = q3 = *sortedBegin;
      iqr = 0.0f;
    }
    else if(N == 2)
    {
      min = *sortedBegin;
      max = *(sortedBegin + 1);
      mean = median = (max + min) * 0.5f;
      q1 = (min + median) * 0.5f;
      q3 = (max + median) * 0.5f;
      iqr = q3 - q1;     
    }
    else
    {
      bool isOddLength = N & 1;
      int half = N * 0.5f;
      int quarter = N * 0.25f;
      int threeFourths = N * 0.75f;
      if (isOddLength)
      {
        q1 = (*(sortedBegin + quarter) + *(sortedBegin + quarter + 1)) * 0.5f;
        median = *(sortedBegin + half);
        q3 = (*(sortedBegin + threeFourths) + *(sortedBegin + threeFourths - 1)) * 0.5f;
      }
      else
      {
        q1 = *(sortedBegin + quarter);
        median = (*(sortedBegin + half) + *(sortedBegin + half - 1)) * 0.5f;
        q3 = *(sortedBegin + threeFourths);
      }
      iqr = q3 - q1;
      min = *sortedBegin;
      max = *(sortedEnd-1);

      mean = 0.0f;
      std::for_each(sortedBegin, sortedEnd, [&mean] (float a) {
        mean += a;
      });
      mean = mean / static_cast<float>(N);
    }

    return std::map<std::string, std::string>
    {
      {LogStrings::EvolutionSectionStrings::mMeanString,    fmt::format("{:g}", mean)   },
      {LogStrings::EvolutionSectionStrings::mMinString,     fmt::format("{:g}", min)    },
      {LogStrings::EvolutionSectionStrings::mMaxString,     fmt::format("{:g}", max)    },
      {LogStrings::EvolutionSectionStrings::mQ1String,      fmt::format("{:g}", q1)     },
      {LogStrings::EvolutionSectionStrings::mMedianString,  fmt::format("{:g}", median) },
      {LogStrings::EvolutionSectionStrings::mQ3String,      fmt::format("{:g}", q3)     },
      {LogStrings::EvolutionSectionStrings::mIqrString,     fmt::format("{:g}", iqr)    }
    };
  }

  template<typename RandomIt>
  std::map<std::string, std::string> getStatistics(const RandomIt sortedBegin,
                                                   const RandomIt sortedEnd) const
  {
    typedef typename std::iterator_traits<RandomIt>::iterator_category category;
    return getStatisticsImpl(sortedBegin, sortedEnd, category());
  }

  std::string serializeVector(const std::vector<float>& vectorToPrint, 
                              std::string separator = LogStrings::CSVStrings::mArraySeparator) const;

  const std::string outputFormatToString() const;
};
#endif //LOGGING_HOOKS_H_INCLUDED