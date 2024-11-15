/**
 * @file        utils.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Source code file containing all definition for auxilliar classes used during the simulation.
 * 
 * @version     0.3
 * 
 * @date        2020-03-02 (created) \n
 *              2020-03-06 (update) \n
 */

#include "utils.h"
#include <iostream>

/**
 * @brief Construct a new auxilliary GridInfo struct
 * 
 * @param [in] xSize    - size of the grid on the x axis
 * @param [in] ySize    - size of the grid on the y axis
 * @param [in] dX       - the distance between the point on the grid in meters
 */
GridInfo::GridInfo(int xSize, 
           int ySize, 
           float dX): 
  nx(xSize), 
  ny(ySize), 
  dx(dX) {}
// END of Grid::ctor

/**
 * @brief Construct a new auxilliary Sonication object
 * 
 * @param [in] xP     - the x coordinate for the sonication
 * @param [in] yP     - the y coordinate for the sonication
 * @param [in] tOn    - how long the sonication runs in seconds
 * @param [in] tOff   - for how long to cool down after the sonication
 * @param [in] dt     - the step size for the simulation
 */
Sonication::Sonication(float xP, 
                       float yP, 
                       float tOn, 
                       float tOff, 
                       float dt):
  xPos(xP),
  yPos(yP),
  timeOn(tOn),
  timeOff(tOff),
  dT(dt) {}
// END of Sonication::ctor

/**
 * @brief Get the discrete amount of steps for how long the sonication runs
 * 
 * @return int    - the number of steps
 */
int Sonication::getOnSteps()
{
  return int(std::round(timeOn/dT));
}
// END of Sonication::getOnSteps

/**
 * @brief Get the discrete amount of steps for how long the medium cools down after sonication
 * 
 * @return int    - the number of steps
 */
int Sonication::getOffSteps()
{
  return int(std::round(timeOff/dT));
}
// END of Sonication::getOffSteps

/**
 * @brief Get the discrete x coordinate for the sonication in the matrix
 * 
 * @return int    - the x coordinate
 */
int Sonication::getXPos()
{
  return int(std::round(xPos));
}
// END of Sonication::getXPos

/**
 * @brief Get the discrete y coordinate for the sonication in the matrix
 * 
 * @return int    - the y coordinate
 */
int Sonication::getYPos()
{
  return int(std::round(yPos));
}
// END of Sonication::getYPos

/**
 * @brief Construct a new Medium object representing the medium on which the sonication is done. \n
 *        Achieved by multiple matrices, each giving a value to a point on the grid.
 * 
 * @param [in] folderPath   - the path to the folder containing density.txt, heat.txt and conductivity.txt
 */
Medium::Medium(std::string folderPath):
  GridInfo(0, 0, 0.0f)
{
  init(folderPath);
}
// END of Medium::ctor

/**
 * @brief Initializes the medium to a begining state. Prepares all data for simulation.
 * 
 * @param [in] folderPath   - the path to the folder containing mediumSpecs.dat, density.dat, heat.dat and conductivity.dat
 */
void Medium::init(std::string folderPath)
{
  if (folderPath.back() != SEP)
  {
    folderPath += SEP;
  }

  std::fstream densityFile(       folderPath + "density.dat",        std::ios_base::in);
  std::fstream specificHeatFile(  folderPath + "heat.dat",           std::ios_base::in);
  std::fstream conductivityFile(  folderPath + "conductivity.dat",   std::ios_base::in);
  std::fstream mediumSpecsFile(   folderPath + "mediumSpecs.dat",    std::ios_base::in);

  if (densityFile.bad() || specificHeatFile.bad() || conductivityFile.bad() || mediumSpecsFile.bad())
  {
    throw std::runtime_error("One or more of the needed medium specifying files could not be opened.");
  }
  
  float              density,          specificHeat,     thermalConductivity;
  std::string        xDimSizeString,   yDimSizeString,   dxSizeString;          // for parsing the mediumSpecs file
  std::string        densityLine,      heatLine,         conductivityLine;      // for parsing lines of data files
  std::string        densityVal,       heatVal,          conductivityVal;       // for single values on each line
  std::stringstream  densityStream,    heatStream,       conductivityStream;    // to convert lines of files into streams

  std::getline(mediumSpecsFile, xDimSizeString, ',');   // nx
  std::getline(mediumSpecsFile, yDimSizeString, ',');   // ny
  std::getline(mediumSpecsFile, dxSizeString, ',');     // dx

  nx = int(strtof(xDimSizeString.c_str(), nullptr));
  ny = int(strtof(yDimSizeString.c_str(), nullptr));
  dx = strtof(dxSizeString.c_str(), nullptr);
  mediumSpecsFile.close();

  // comma operator to perform every reading and avoid shortcircuiting on first
  while(std::getline(densityFile,      densityLine),
        std::getline(specificHeatFile, heatLine),
        std::getline(conductivityFile, conductivityLine))
  {
    densityStream.clear();
    heatStream.clear();
    conductivityStream.clear();

    densityStream.str(densityLine);
    heatStream.str(heatLine);
    conductivityStream.str(conductivityLine);

    densityVal      = "";
    heatVal         = "";
    conductivityVal = "";

    while(std::getline(densityStream,      densityVal,      ','),
          std::getline(heatStream,         heatVal,         ','),
          std::getline(conductivityStream, conductivityVal, ','))
    {
      density             = strtof(densityVal.c_str(), nullptr);
      specificHeat        = strtof(heatVal.c_str(), nullptr);
      thermalConductivity = strtof(conductivityVal.c_str(), nullptr);
      
      q.push_back(0.0f);                          // volume rate of heat deposition [W/m^3]
      t.push_back(37.0f);                         // temperature distribution [degC]
      rho.push_back(density);                     // tissue mass density [kg/m^3]
      c.push_back(specificHeat);                  // tissue specific heat capacity [J/(kg.K)]
      cem43.push_back(0.0f);                      // thermal dose given in cumulative equivalent minutes (cem) relative to T = 43 degC
      lambda.push_back(thermalConductivity);      // tissue thermal conductivity [W/(m.K)]
    }
  }

  if (!(densityFile.eof() && specificHeatFile.eof() && conductivityFile.eof()))
  {
    throw std::runtime_error("One or more of the needed medium files are not equal in their dimensions.");
  }
  if (densityFile.bad() || specificHeatFile.bad() || conductivityFile.bad())
  {
    throw std::runtime_error("Error reading medium defining files.");
  }
}
// END of Medium::init