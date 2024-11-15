/**
 * @file        calculateHeating.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file for creating target maps describing the environment
 * 
 * @version     0.1
 * 
 * @date        2020-03-20 (created) \n
 */
#pragma once
#ifndef MAPS_FACTORY_H
#define MAPS_FACTORY_H

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

#include "discreteMap.h"

/**
 * @brief Creates a lesion map from cem43 map
 * 
 * @param [in] cem43          - a stl vector container with cem43 data from heating equations
 * @return std::vector<int>   - a binary lesion map (1 if cem43 .> 240 else 0)
 */
std::vector<int> CreateLesionMap(std::vector<int>& cem43);
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
                          std::function<int(int, int)> ftor
                          );
#endif