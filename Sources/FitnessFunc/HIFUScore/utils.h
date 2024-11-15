/**
 * @file        utils.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Header file containing all declarations for auxilliar classes used during the simulation.
 * 
 * @version     0.3
 * 
 * @date        2020-03-02 (created) \n
 *              2020-03-06 (update) \n
 */
#pragma once
#ifndef UTILS_STRUCTS_H
#define UTILS_STRUCTS_H

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#ifndef SEP
  #ifdef _WIN32
    #define SEP '\\';
  #else 
    #define SEP '/'
  #endif
#endif

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
  unsigned nx;

  /**
   * @brief Size of the grid on the y axis
   * 
   */
  unsigned ny;

  /**
   * @brief The distance between the point on the grid in meters
   * 
   */
  float dx;
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
           float dX);

};
//END of Grid struct

/**
 * @brief An auxilliary sonication object carrying all the neccessary data \n
 *        to start a sonication simulation given a medium.
 * 
 */
struct Sonication
{
  /**
   * @brief The x coordinate for the sonication
   * 
   */
  const float xPos;

  /**
   * @brief The y coordinate for the sonication
   * 
   */
  const float yPos;

  /**
   * @brief For how long the sonication runs
   * 
   */
  const float timeOn;

  /**
   * @brief For how long to cool down after the sonication
   * 
   */
  const float timeOff;

  /**
   * @brief The step size for the simulation
   * 
   */
  const float dT;

  /**
   * @brief Construct a new auxilliary Sonication object
   * 
   * @param [in] xP     - the x coordinate for the sonication
   * @param [in] yP     - the y coordinate for the sonication
   * @param [in] tOn    - for how long the sonication runs
   * @param [in] tOff   - for how long to cool down after the sonication
   * @param [in] dt     - the step size for the simulation
   */
  Sonication(float xP,
             float yP,
             float tOn,
             float tOff,
             float dt);

  /**
   * @brief Get the discrete amount of steps for how long the sonication runs
   * 
   * @return int    - the number of steps
   */            
  int getOnSteps();

  /**
   * @brief Get the discrete amount of steps for how long the medium cools down after sonication
   * 
   * @return int    - the number of steps
   */
  int getOffSteps();

  /**
   * @brief Get the discrete x coordinate for the sonication in the matrix
   * 
   * @return int    - the x coordinate
   */
  int getXPos();

  /**
   * @brief Get the discrete y coordinate for the sonication in the matrix
   * 
   * @return int    - the y coordinate
   */
  int getYPos();
};
//END of Sonication struct

/**
 * @brief A class representing the medium on which the heat diffusion simulation is performed \n
 *        Inherits from Grid class that represents the informations about the size of the \n
 *        discrete grid.
 */
class Medium : public GridInfo
{
private:
  /**
   * @brief Initializes the medium to a begining state. Prepares all data for simulation.
   * 
   * @param [in] folderPath   - the path to the folder containing density.txt, heat.txt and conductivity.txt
   */
  void init(std::string folderPath);

public:
  /**
   * @brief volume rate of heat deposition [W/m^3]
   * 
   */
  std::vector<float> q;

  /**
   * @brief temperature distribution [degC]
   * 
   */
  std::vector<float> t;

  /**
   * @brief tissue mass density [kg/m^3]
   * 
   */
  std::vector<float> rho;

  /**
   * @brief tissue specific heat capacity [J/(kg.K)]
   * 
   */
  std::vector<float> c;

  /**
   * @brief thermal dose given in cumulative equivalent minutes (cem) relative to T = 43 degC
   * 
   */
  std::vector<float> cem43;

  /**
   * @brief tissue thermal conductivity [W/(m.K)]
   * 
   */
  std::vector<float> lambda;
  
  /**
   * @brief Construct a new Medium object representing the medium on which the sonication is done. \n
   *        Achieved by multiple matrices, each giving a value to a point on the grid.
   * 
   * @param [in] folderPath   - the path to the folder containing density.txt, heat.txt and conductivity.txt
   */
  Medium(std::string folderPath);

};
//END of Medium class
#endif