/**
 * @file        calculateHeating.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for creating target maps describing the environment
 * 
 * @version     0.1
 * 
 * @date        2020-03-20 (created) \n
 */

#include "mapsFactory.h"

/**
 * @brief Creates a lesion map from cem43 map
 * 
 * @param [in] cem43          - a stl vector container with cem43 data from heating equations
 * @return std::vector<int>   - a binary lesion map (1 if cem43 .>= 240 else 0)
 */
std::vector<int> CreateLesionMap(std::vector<int>& cem43)
{
  std::vector<int> lesionMap;
  auto filter = [&lesionMap](int val) { (val >= 240) ? lesionMap.push_back(1) : lesionMap.push_back(0); };
  std::for_each(cem43.begin(), cem43.end(), filter);
  return lesionMap;
}
//END of CreateLesionMap

/**
 * @brief Creates a DiscreteMap object of the required size from the circles defining file
 * 
 * @param [in] nx                     - size of the map on x dimension
 * @param [in] ny                     - size of the map on y dimension
 * @param [in] dx                     - distance between the neighbours in the real medium [m]
 * @param [in] circlesMapFilePath     - path to the file containing the definitions of circles that creates the \n
 *                                      discrete map
 * @return DiscreteMap                - a nx x ny map filled with circles given in the file
 */
DiscreteMap GetCirclesMap(const unsigned nx,
                          const unsigned ny,
                          const float dx,
                          std::string circlesMapFilePath,
                          std::function<int(int, int)> ftor)
{
  std::vector<CircleInfo> circlesCollection;
  std::string circleLine;
  std::stringstream lineStream;

  std::fstream mapCirclesStream(circlesMapFilePath, std::ios_base::in);
  if (mapCirclesStream.bad())
  {
    throw std::runtime_error("File " + circlesMapFilePath + " with circles specification does not exist");
  }

  std::vector<DiscreteMap> maps;

  std::string sx, sy, radius, val;

  try
  {
    while (std::getline(mapCirclesStream, circleLine))
    {
      lineStream.str(circleLine);
      if (lineStream.peek() != '#')
      {
        if (lineStream.peek() == '-')
        {
          maps.push_back(DiscreteMap::CreateTargetsMap(nx, ny, dx, circlesCollection, ftor));
          circlesCollection.clear();
        }
        else
        {
          std::getline(lineStream,  sx,      ',');
          std::getline(lineStream,  sy,      ',');
          std::getline(lineStream,  radius,  ',');
          std::getline(lineStream,  val,     ',');

          auto n = val.find("#");
          if (n != std::string::npos)
          {
            val = val.substr(0, n);
          }

          circlesCollection.push_back(
            CircleInfo { 
              unsigned(std::round(strtof(sx.c_str(), nullptr))),
              unsigned(std::round(strtof(sy.c_str(), nullptr))),
              unsigned(std::round(strtof(radius.c_str(), nullptr))),
              int(std::round(strtof(val.c_str(), nullptr))),
            }
          );
        }
      }
      lineStream.clear();
    }
  }
  catch(...)
  {
    throw std::runtime_error("File " + circlesMapFilePath + " contains invalid characters.");
  }
  
  if (!circlesCollection.empty())
  {
    maps.push_back(DiscreteMap::CreateTargetsMap(nx, ny, dx, circlesCollection, ftor));
    circlesCollection.clear();
  }

  for (int i = maps.size()-1; i > 0; i--)
  {
    maps[i-1] = maps[i] + maps[i-1];
  }
  return maps[0];
}
//END of GetCirclesMap