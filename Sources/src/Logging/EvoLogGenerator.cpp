/**
 * @file        EvoLogGenerator.cpp
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

#include <Logging/EvoLogGenerator.hpp>
#include <Logging/LogStrings.hpp>

#include <fmt/format.h>
#include <nlohmann/json.hpp>

#include <vector>
#include <map>
#include <deque>
#include <iterator>
#include <string>
#include <algorithm>  //find_if
#include <chrono>
#include <stdexcept>


void CustomLogData::add(const std::string& key, 
                        const std::string& val)
{
  mCustomData.push_back(KeyValuePair(key, val));
}

void CustomLogData::add(const std::pair<std::string, std::string>& stlPair)
{
  mCustomData.push_back(KeyValuePair(stlPair));
}

KeyValuePair& CustomLogData::operator[](const int& index)
{
  return mCustomData[index];
}

auto CustomLogData::begin()         { return mCustomData.begin();  }
auto CustomLogData::end()           { return mCustomData.end();    }
auto CustomLogData::begin()  const  { return mCustomData.begin();  }
auto CustomLogData::end()    const  { return mCustomData.end();    }
auto CustomLogData::cbegin() const  { return mCustomData.cbegin(); }
auto CustomLogData::cend()   const  { return mCustomData.cend();   }
auto CustomLogData::size()   const  { return mCustomData.size();   }
auto CustomLogData::empty()  const  { return mCustomData.empty();  }
  
CustomLogData::CustomLogData (const std::map<std::string, std::string>& stlMap)
{
  for (auto& dataPair : stlMap)
  {
    mCustomData.push_back(std::pair<std::string, std::string>(dataPair.first, dataPair.second));
  }
}

CustomLogData::CustomLogData (const std::vector<std::pair<std::string, std::string>>& stlVector)
{
  for (auto& dataPair : stlVector)
  {
    mCustomData.push_back(std::pair<std::string, std::string>(dataPair.first, dataPair.second));
  }
}

CustomLogData::CustomLogData(const CustomLogData& other) :
  mCustomData(other.mCustomData)
{}

CustomLogData::CustomLogData()  {}



EvoLogGenerator::EvoLogGenerator(OutputFormat outputFormat):
  mOutputFormat(outputFormat)
{}

std::string 
EvoLogGenerator::createLaunchLogString(const int&                 popSize,
                                       const int&                 dimension,
                                       const long long&           seed,
                                       const std::time_t          dateTime,
                                       const std::vector<float>&  upperBounds,
                                       const std::vector<float>&  lowerBounds,                                    
                                       const CustomLogData&       customData)
{
  auto dataMapForOrdering = getLaunchMap(popSize, dimension, seed, dateTime, upperBounds, lowerBounds);
  auto pairsInOrder       = orderRequiredDataAndAddCustom(dataMapForOrdering, LogStrings::LaunchSectionStrings::mOrder, customData);

  std::string retVal("");
  switch (mOutputFormat)
  {
    case OutputFormat::CSV:
      retVal = createCSVString(pairsInOrder, true);
      break;
    case OutputFormat::JSON:
      retVal = createJSONString(pairsInOrder);
      break;
    case OutputFormat::SimpleParse:
      retVal = createSimpleParseString(pairsInOrder);
      break;
  }
  return retVal;
}

std::string 
EvoLogGenerator::createEvolutionLogString(const int&                currentGeneration,
                                          const std::vector<float>& populationFitnesses,
                                          const float&              bestSoFar,
                                          const CustomLogData&      customData)
{

  auto dataMapForOrdering = getGenerationMap(currentGeneration, populationFitnesses, -0.0f, bestSoFar, false);
  auto pairsInOrder       = orderRequiredDataAndAddCustom(dataMapForOrdering, LogStrings::EvolutionSectionStrings::mOrder, customData);

  std::string retVal("");
  switch (mOutputFormat)
  {
    case OutputFormat::CSV:
      retVal = createCSVString(pairsInOrder, mCreateCSVGenHeader);
      mCreateCSVGenHeader = false;
      break;
    case OutputFormat::JSON:
      retVal = createJSONString(pairsInOrder);
      break;
    case OutputFormat::SimpleParse:
      retVal = createSimpleParseString(pairsInOrder);
      break;
  }
  return retVal;
}


std::string 
EvoLogGenerator::createEvolutionLogString(const int&            currentGeneration,
                                          const float&          currentFitness,
                                          const float&          bestSoFar,
                                          const CustomLogData&  customData)
{
  auto dataMapForOrdering = getGenerationMap(currentGeneration, std::vector<float>(), currentFitness, bestSoFar, true);
  auto pairsInOrder       = orderRequiredDataAndAddCustom(dataMapForOrdering, LogStrings::EvolutionSectionStrings::mOrder, customData);

  std::string retVal("");
  switch (mOutputFormat)
  {
    case OutputFormat::CSV:
      retVal = createCSVString(pairsInOrder, mCreateCSVGenHeader);
      mCreateCSVGenHeader = false;
      break;
    case OutputFormat::JSON:
      retVal = createJSONString(pairsInOrder);
      break;
    case OutputFormat::SimpleParse:
      retVal = createSimpleParseString(pairsInOrder);
      break;
  }
  return retVal;
}





std::string 
EvoLogGenerator::createResultsLogString(const float&              bestFitness,
                                        const std::time_t&        runTime,
                                        const int&                fitnessEvaluations,
                                        const int&                generations,
                                        const std::vector<float>& bestChromo,
                                        const CustomLogData&      customData)
{
  auto dataMapForOrdering = getResultsMap(bestFitness, runTime, fitnessEvaluations, generations, bestChromo);
  auto pairsInOrder       = orderRequiredDataAndAddCustom(dataMapForOrdering, LogStrings::ResultSectionStrings::mOrder, customData);

  std::string retVal("");
  switch (mOutputFormat)
  {
    case OutputFormat::CSV:
      retVal = createCSVString(pairsInOrder, true);
      break;
    case OutputFormat::JSON:
      retVal = createJSONString(pairsInOrder);
      break;
    case OutputFormat::SimpleParse:
      retVal = createSimpleParseString(pairsInOrder);
      break;
  }
  return retVal;
}

std::string
EvoLogGenerator::createCSVString(const std::vector<KeyValuePair>& pairsToPrint,
                                 const bool                       printHeader) const
{
  std::string descriptors("");
  std::string values("");

  if (printHeader)
  {
    for (auto& descriptorValuePair : pairsToPrint)
    {
      descriptors += fmt::format("{:s}{:s}", descriptorValuePair.first,  
                                            LogStrings::CSVStrings::mCsvSeparator);
    }
    descriptors += "\n";
  }

  for (auto& descriptorValuePair : pairsToPrint)
  {
    values      += fmt::format("{:s}{:s}", descriptorValuePair.second, 
                                            LogStrings::CSVStrings::mCsvSeparator);
  }
  return fmt::format("{:s}{:s}{:s}", descriptors, values , "\n");
}

std::string 
EvoLogGenerator::createJSONString(const std::vector<KeyValuePair>& pairsToPrint) const
{
  nlohmann::ordered_json header;

  for (auto& descriptorValuePair : pairsToPrint)
  {
    header.emplace(descriptorValuePair.first, descriptorValuePair.second);
  }
  return header.dump(2);
}

std::string 
EvoLogGenerator::createSimpleParseString(const std::vector<KeyValuePair>& pairsToPrint) const
{
  std::string toPrint("");
  for (auto& descriptorValuePair : pairsToPrint)
  {
    toPrint += fmt::format("{:s}{:s}{:s}{:s}",
                            LogStrings::SimpleParseStrings::mDataSeparator,
                            descriptorValuePair.first,
                            LogStrings::SimpleParseStrings::mDescriptorValSeparator,
                            descriptorValuePair.second);
  }
  return fmt::format("{:s}{:s}", toPrint, LogStrings::SimpleParseStrings::mGenSeparator);
}

const std::map<std::string, std::string>  
EvoLogGenerator::getLaunchMap(const int& popSize,
                              const int& dimension,
                              const long long& seed,
                              const std::time_t dateTime,
                              const std::vector<float>& upperBounds,
                              const std::vector<float>& lowerBounds) const
{
  auto tm = std::gmtime(&dateTime);
  auto timeString = fmt::format("{:d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}",  tm->tm_year + 1900, 
                                                                            tm->tm_mon + 1, 
                                                                            tm->tm_mday, 
                                                                            tm->tm_hour,
                                                                            tm->tm_min,
                                                                            tm->tm_sec);
  return std::map<std::string, std::string>
  {
    {LogStrings::LaunchSectionStrings::mDateTimeString,     timeString},
    {LogStrings::LaunchSectionStrings::mPopSizeString,      fmt::format("{:d}", popSize)},
    {LogStrings::LaunchSectionStrings::mDimensionString,    fmt::format("{:d}", dimension)},
    {LogStrings::LaunchSectionStrings::mSeedString,         fmt::format("{:d}", seed)},
    {LogStrings::LaunchSectionStrings::mUpperBoundsString,  fmt::format("{:s}", serializeVector(upperBounds))},
    {LogStrings::LaunchSectionStrings::mLowerBoundsString,  fmt::format("{:s}", serializeVector(lowerBounds))},
    {LogStrings::LaunchSectionStrings::mLogModeString,      std::string(outputFormatToString())} 
  };
}

const std::map<std::string, std::string> 
EvoLogGenerator::getGenerationMap(const int&                currentGeneration,
                                  const std::vector<float>& populationFitnesses,
                                  const float&              currentFitness,
                                  const float&              bestSoFar,
                                  const bool                rolling)
{

  std::map<std::string, std::string> descriptorToValueMap
  {
    {LogStrings::EvolutionSectionStrings::mGenString,      fmt::format("{:d}",  currentGeneration)},
    {LogStrings::EvolutionSectionStrings::mBestString,     fmt::format("{:g}",  bestSoFar)},
    {LogStrings::EvolutionSectionStrings::mPopDumpString,  fmt::format("{:s}",  serializeVector(populationFitnesses))}
  };

  std::map<std::string, std::string> stats;
  if (rolling)
  {
    auto insertIt = std::upper_bound(mHistory.begin(), mHistory.end(), currentFitness);
    mHistory.insert(insertIt, currentFitness);
    stats = getStatistics(mHistory.cbegin(), mHistory.cend());
  }
  else
  {
    std::vector<float> popCopy(populationFitnesses.size(), 0.0f);
    std::copy(populationFitnesses.cbegin(), populationFitnesses.cend(), popCopy.begin());
    std::sort(popCopy.begin(), popCopy.end(), std::less<float>());
    stats = getStatistics(popCopy.cbegin(), popCopy.cend());
  }

  descriptorToValueMap.insert(stats.begin(), stats.end());
  return descriptorToValueMap;
}

const std::map<std::string, std::string> 
EvoLogGenerator::getResultsMap(const float&               bestFitness,
                               const std::time_t&         runTime,
                               const int&                 fitnessEvaluations,
                               const int&                 generations,
                               const std::vector<float>&  bestResult) const
{
  return std::map<std::string, std::string>
  {
    {LogStrings::ResultSectionStrings::mResFitnessString,   fmt::format("{:g}", bestFitness)},
    {LogStrings::ResultSectionStrings::mTimeElapsedString,  fmt::format("{:d}", static_cast<long long>(runTime))},
    {LogStrings::ResultSectionStrings::mFitEvalstring,      fmt::format("{:d}", fitnessEvaluations)},
    {LogStrings::ResultSectionStrings::mGensNeededString,   fmt::format("{:d}", generations)},
    {LogStrings::ResultSectionStrings::mBestString,         fmt::format("{:s}", serializeVector(bestResult))}
  };
}

const std::vector<KeyValuePair> 
EvoLogGenerator::orderRequiredDataAndAddCustom(const std::map<std::string, std::string> requiredData,
                                               const std::vector<std::string>&          orderingOfRequiredData,
                                               const CustomLogData&                     customData) const
{
  auto data = std::vector<KeyValuePair>();
  for (auto& key : orderingOfRequiredData)
  {
    if (requiredData.find(key) != requiredData.cend())
    {
      data.push_back({key, requiredData.at(key)});
    }
    else
    {
      data.push_back({key, "UNSPECIFIED"});
    }
  }

  for (auto& userPair: customData)
  {
    auto pairIndex = std::find_if(data.begin(), data.end(), 
      [&userPair] (const KeyValuePair& val) -> bool
      {
        return val.first == userPair.first;
      }
    );

    if (pairIndex != data.end())
    {
      data[std::distance(data.begin(),pairIndex)].first = userPair.first;
    }
    else
    {
      data.push_back(userPair);
    }
  }
  return data;
}

std::string 
EvoLogGenerator::serializeVector(const  std::vector<float>&  vectorToPrint, 
                                        std::string          separator) const
{
  return fmt::format("{:g}", fmt::join(vectorToPrint, separator));
}

const std::string 
EvoLogGenerator::outputFormatToString() const
{
  switch (mOutputFormat)
  {
    case OutputFormat::CSV:
      return "CSV";
    case OutputFormat::JSON:
      return "JSON";
    case OutputFormat::SimpleParse:
      return "SimpleParse";
    default:
      return "UNSPECIFIED";
  };
}
