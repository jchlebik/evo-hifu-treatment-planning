#pragma once
#ifndef KEVOOPT_FITNESSFUNCTIONS_HIFUSCORE_DATAMODELS_GRIDINFO_HPP_INCLUDE
#define KEVOOPT_FITNESSFUNCTIONS_HIFUSCORE_DATAMODELS_GRIDINFO_HPP_INCLUDE

/**
 * @brief A structure carrying information about the discrete grid \n
 *        representing the search space
 * 
 */
class GridInfo
{
public:
  /**
   * @brief Size of the grid on the x axis
   * 
   */
  unsigned mNx;

  /**
   * @brief Size of the grid on the y axis
   * 
   */
  unsigned mNy;

  /**
   * @brief The distance between the point on the grid in meters
   * 
   */
  float mDx;
protected:
  /**
   * @brief Construct a new auxilliary Grid object
   * 
   * @param [in] xSize    - size of the grid on the x axis
   * @param [in] ySize    - size of the grid on the y axis
   * @param [in] dX       - the distance between the point on the grid in meters
   */
  GridInfo(int xSize,
           int ySize,
           float dX):
  mNx(xSize), 
  mNy(ySize), 
  mDx(dX) {}
// END of Grid::ctor

};
//END of Grid struct
#endif