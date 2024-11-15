/**
 * @file        Medium.cpp
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Source file containing all definitions for Medium class used during the simulation. \n
 *              This class represents and contains all the data about the medium we are using for \n
 *              heat deposit and tissue ablation simulations.
 * 
 * @version     0.1
 * 
 * @date        2020-11-24 (created) \n
 */

#include <FitnessFunctions/HIFUScore/DataModels/Medium.hpp>

#include <FitnessFunctions/HIFUScore/DataModels/utils.hpp>
#include <FitnessFunctions/HIFUScore/DataModels/GridInfo.hpp>
#include <FitnessFunctions/HIFUScore/DataModels/HeatTargetWindow.hpp>

#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <string>

#include <cstdlib>

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

Medium::~Medium()
{
  _mm_free(mQ);
  _mm_free(mT);
  _mm_free(mCem43);
  _mm_free(mRho);
  _mm_free(mC);
  _mm_free(mLambda);
}

/**
 * @brief Construct a new Medium object representing the medium on which the sonication is done. \n
 *        Achieved by multiple matrices, each giving a value to a point on the grid and cropped by \n
 *        a given window represented by HeatTargetWindow class.
 * 
 * @param [in] targetWindow   - The instance of the HeatTargetWindow
 */
Medium::Medium(HeatTargetWindow   targetWindow, 
               float* heat, 
               float* temperature, 
               float* density, 
               float* capacity,
               float* cem,
               float* conductivity) :
  GridInfo(targetWindow.getXAxisLen(), targetWindow.getYAxisLen(), targetWindow.getResolution())
{
  mMatSize = mNx * mNy;
  mQ       = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));
  mT       = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));
  mCem43   = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));
  mRho     = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));  
  mC       = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ )); 
  mLambda  = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ )); 


  for (int i = 0; i < mMatSize; i++)
  {
    mQ[i] = heat[i];
    mT[i] = temperature[i];
    mCem43[i] = cem[i];
    mRho[i] = density[i];
    mC[i] = capacity[i];
    mLambda[i] = conductivity[i];
  }
}
//END of Medium::ctor


Medium::Medium(const Medium & cpy) :
  GridInfo(cpy.mNx, cpy.mNy, cpy.mDx)
{
  mMatSize = mNx * mNy;
  mQ       = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));
  mT       = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));
  mCem43   = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));
  mRho     = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));
  mC       = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ )); 
  mLambda  = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ )); 

  for (int i = 0; i < mMatSize; i++)
  {
    mQ[i]       = cpy.mQ[i];
    mT[i]       = cpy.mT[i];
    mCem43[i]   = cpy.mCem43[i];
    mRho[i]     = cpy.mRho[i];
    mC[i]       = cpy.mC[i];
    mLambda[i]  = cpy.mLambda[i];
  }
}


/**
 * @brief Initializes the medium to a begining state. Prepares all data for simulation.
 * 
 * @param [in] folderPath   - the path to the folder containing mediumSpecs.dat, density.dat, heat.dat and conductivity.dat
 */
void Medium::init(std::string folderPath)
{
  if ( folderPath.back() != SEP)
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
  std::getline(mediumSpecsFile, dxSizeString,   ',');   // dx

  mNx = int(std::strtof(xDimSizeString.c_str(), nullptr));
  mNy = int(std::strtof(yDimSizeString.c_str(), nullptr));
  mDx = std::strtof(dxSizeString.c_str(), nullptr);
  mediumSpecsFile.close();

  mMatSize = mNx * mNy;
  mQ       = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));
  mT       = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));
  mCem43   = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));
  mRho     = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ ));  
  mC       = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ )); 
  mLambda  = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * mMatSize, _ALIGN_ )); 

  int i = 0;
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
      
      mQ[i]      = 0.0f;                          // volume rate of heat deposition [W/m^3]
      mT[i]      = 37.0f;                         // temperature distribution [degC]
      mCem43[i]  = 0.0f;                      // thermal dose given in cumulative equivalent minutes (cem) relative to T = 43 degC
      mRho[i]    = density;                     // tissue mass density [kg/m^3]
      mC[i]      = specificHeat;                  // tissue specific heat capacity [J/(kg.K)]
      mLambda[i] = thermalConductivity;      // tissue thermal conductivity [W/(m.K)]
      i++;
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



void Medium::setHeatDeposit(float* heatAdded)
{
  HeatTargetWindow wholeMedium = HeatTargetWindow(0, 0, mNx, mNy, mNx, mNy, mDx);
  setHeatDepositEnergyOnWindow(wholeMedium, heatAdded);
}

/**
 * @brief Sets the heat deposit energy into the medium after cropping it by the provided window.
 * 
 * @param [in] heatTargetWindow       - Instance of HeatTargetWindow used to crop the heat doposition \n
 *                                      vector to the appropriate window.
 * @param [in] heatAdded              - Vector of data containing the values for the added heat energy \n
 *                                      onto the medium. Expected to be of the original size and to be \n
 *                                      cropped by the provided window.
 */
void Medium::setHeatDepositEnergyOnWindow(HeatTargetWindow & heatTargetWindow, float* heatAdded)
{
  const unsigned xRange = heatTargetWindow.getXAxisLen();
  const unsigned yRange = heatTargetWindow.getYAxisLen();
  const unsigned xStart = heatTargetWindow.getXStart();
  const unsigned yStart = heatTargetWindow.getYStart();

  const unsigned xRangeOriginal = heatTargetWindow.getXAxisOriginalLen();

  for (int row = 0; row < yRange; row++)
  {
    const int rowStartIndex = (yStart + row) * xRangeOriginal + xStart;
    for(int col = 0; col < xRange; col++)
    {
      mQ[row * xRange + col] = heatAdded[rowStartIndex + col];
    }
  }
}
//END of setHeatDepositEnergyOnWindow

/**
 * @brief Static method to create an original sized medium with updated values \n
 *        at locations in the medium cropped by window. Does not modify its \n
 *        parameters, instead opting to create a new instance.
 * 
 * @param [in] heatTargetWindow   - A window that was used to create the croppedMedium \n
 *                                  instance.
 * @param [in] originalMedium     - Original sized Medium containing the data from before \n 
 *                                  the sonications. Copy of this medium with updated values \n
 *                                  for t and cem43 vector is returned.
 * @param [in] croppedMedium      - A medium that was cropped by the provided window carrying \n
 *                                  all relevant information after the sonications.
 * 
 * @returns Medium                - A copy of the originalMedium parameter with updated values \n
 *                                  of t and cem43 vectors from the cropped medium.
 * 
 */
Medium Medium::updateOriginalByCroppedMedium(HeatTargetWindow & heatTargetWindow,
                                             Medium originalMedium,
                                             Medium & croppedMedium)
{
  const int xRange = heatTargetWindow.getXAxisLen();
  const int yRange = heatTargetWindow.getYAxisLen();
  const int xStart = heatTargetWindow.getXStart();
  const int yStart = heatTargetWindow.getYStart();

  Medium md(originalMedium);

  for (int row = 0; row < yRange; row++)
  {
    const int rowStartIndex = (yStart + row) * originalMedium.mNx + xStart;
    for (int col = 0; col < xRange; col++)
    {
      md.mT[rowStartIndex + col] = croppedMedium.mT[row * xRange + col];
      md.mCem43[rowStartIndex + col] = croppedMedium.mCem43[row * xRange + col];
    }
  }

  return md;
}
//END of updateOriginalByCroppedMedium
/**
 * @brief   Provides a new instance of the Medium class that contains the data \n
 *          that fit inside the provided HeatTargetWindow. The new instance is \n
 *          of smaller size (the size of the window).
 * 
 * @param [in] heatTargetWindow     - A window with which to crop the medium.
 * 
 * @returns Medium                  - A new smaller instance of the Medium class \n
 *                                    with data only inside the window.
 * 
 */
Medium Medium::getWindowedMedium(HeatTargetWindow & heatTargetWindow)
{
  const int xRange = heatTargetWindow.getXAxisLen();
  const int yRange = heatTargetWindow.getYAxisLen();
  const int xStart = heatTargetWindow.getXStart();
  const int yStart = heatTargetWindow.getYStart();

  const unsigned xRangeOriginal = heatTargetWindow.getXAxisOriginalLen();

  std::vector<float> heat(xRange * yRange);
  std::vector<float> temperature(xRange * yRange);
  std::vector<float> density(xRange * yRange);
  std::vector<float> capacity(xRange * yRange);
  std::vector<float> cem(xRange * yRange);
  std::vector<float> conductivity(xRange * yRange);
  
  
  std::vector<std::vector<float>*> croppedDatasVector
  {
    &heat, &temperature, &density, &capacity, &cem, &conductivity
  };

  std::vector<float**> originalDatasVector
  {
    &mQ, &mT, &mRho, &mC, &mCem43, &mLambda
  };

  for (int i = 0; i < originalDatasVector.size(); i++)
  {
    for (int row = 0; row < yRange; row++)
    {
      const int rowStartIndex = (yStart + row) * xRangeOriginal + xStart;
      for(int col = 0; col < xRange; col++)
      {
        (*croppedDatasVector[i])[row * xRange + col] = (*originalDatasVector[i])[rowStartIndex + col];
      }
    }
  }

  return Medium(heatTargetWindow,
                heat.data(),
                temperature.data(),
                density.data(),
                capacity.data(),
                cem.data(),
                conductivity.data());
}
//END of getWindowedMedium
