/**
 * @file        discreteMap.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Header file containing all declarations for DiscreteMap class and its data structure.
 * 
 * @version     0.2
 * 
 * @date        2020-03-12 (created) \n
 *              2020-03-18 (update) \n
 */
#pragma once
#ifndef TARGET_MAP_CLASS_H
#define TARGET_MAP_CLASS_H

#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <stdlib.h>
#include "utils.h"


/**
 * @brief An auxilliary data structure to wrap the definition of a circle in a target map
 * 
 */
struct CircleInfo
{
  /**
   * @brief an x coordinate in the grid of the circle center
   * 
   */
  const unsigned x;
  /**
   * @brief an y coordinate in the grid of the circle center
   * 
   */
  const unsigned y;
  /**
   * @brief a radius of the circle
   * 
   */
  const unsigned radius;
  /**
   * @brief a value that the circle is going to be filled with
   * 
   */
  const int      circleValue;
};

/**
 * @brief A class representing the discrete map of values. 
 * 
 */
class DiscreteMap : public GridInfo
{
  /**
   * @brief a vector of values containing the disrete map of values
   * 
   */
  std::vector<int> data;

public:
 /**
  * @brief Construct a new 2D TargetMap object of given Grid size with circles specified in coordinates. \n
  *        Edge cases are wrapped around its opposite edge calculated by modulus operation. Identical to the 
  *        Matlab version used for fitness maps in k-Wave EUD.
  * 
  * @param [in] g              - a GridInfo struct carrying relevant information about the map (xn, yn, dx)
  * @param [in] data           - the flattened data values of the map.
  */  
  DiscreteMap(GridInfo g, 
              std::vector<int>& mapData);

  /**
   * @brief Construct a new 2D TargetMap object of given Grid size with circles specified in coordinates. \n
   *        Edge cases are wrapped around its opposite edge calculated by modulus operation. Identical to the 
   *        Matlab version used for fitness maps in k-Wave EUD.
   * 
   * @param [in] xDimSize       - the size of the x axis of the discrete medium space
   * @param [in] yDimSize       - the size of the y axis of the discrete medium space
   * @param [in] dx             - the distance between two points in the discrete space in meters
   * @param [in] data           - the flattened data values of the map.
   */  
  DiscreteMap(const unsigned xDimSize,
              const unsigned yDimSize,
              const float dx,
              std::vector<int>& mapData);

  /**
   * @brief sums the data contained in the map
   * 
   * @return int - sum of the data in the map
   */
  int sum();

  /**
   * @brief Creates a lesion map specific to the current state of the given medium
   * 
   * @param [in] m             - Medium instance to create the lesion map for
   * @return DiscreteMap       - a binary DiscreteMap representing the lesion map
   */
  static DiscreteMap CreateLesionMap(Medium& m);

  /**
   * @brief Creates a DiscreteMap object of the required size from the circles defining file
   * 
   * @param [in] nx                     - size of the map on x dimension
   * @param [in] ny                     - size of the map on y dimension
   * @param [in] dx                     - distance between the neighbours in the real medium [m]
   * @param [in] circles                - a vector container of the CircleInfo instances representing the circles on the \n
   *                                      discrete map
   * @return DiscreteMap                - a nx x ny map filled with circles given in container
   */
  static DiscreteMap CreateTargetsMap(const unsigned nx, 
                                      const unsigned ny, 
                                      const float dx, 
                                      std::vector<CircleInfo>& circles,
                                      std::function<int(int, int)> ftor);

  void PrintMap(std::string sep = " ")
  {
    for (int i = 0; i < ny; i++)
    {
      for (int j = 0; j < nx; j++)
      {
        std::cout<<data[i*ny + j] << sep;
      }
      std::cout << std::endl;
    }
  }

  /**
   * @brief Definition of the substration operation for the discrete maps. \n
   *        Equivalent to substraction of matrices. Requires L an R to \n
   *        be of the same size
   * 
   * @param [in] L          - L side argument of the substraction
   * @param [in] R          - R side argument of the substraction
   * @return DiscreteMap    - L - R
   */
  friend DiscreteMap operator-(DiscreteMap const &, DiscreteMap const &);

  /**
   * @brief Definition of the multiplication operation of the discrete maps. \n
   *        Equivalent to scalar product of matrices. Requires L an R to \n
   *        be of the same size
   * 
   * @param [in] L          - L side argument of the multiplication
   * @param [in] R          - R side argument of the multiplication
   * @return DiscreteMap    - L * R
   */
  friend DiscreteMap operator*(DiscreteMap const &, DiscreteMap const &);

  /**
   * @brief Definition of the addition operation of the discrete maps. \n
   *        Equivalent to dot addition of matrices. Requires L an R to \n
   *        be of the same size
   * 
   * @param [in] L          - L side argument of the addition
   * @param [in] R          - R side argument of the addition
   * @return DiscreteMap    - L .+ R
   */
  friend DiscreteMap operator+(DiscreteMap const &, DiscreteMap const &);
};
#endif