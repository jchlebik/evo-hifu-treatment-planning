/**
 * @file        Sonication.hpp
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
 *              2020-11-24 (update) \n
 */

#pragma once
#ifndef KEVOOPT_FITNESSFUNCTIONS_HIFUSCORE_DATAMODELS_SONICATION_HPP_INCLUDE
#define KEVOOPT_FITNESSFUNCTIONS_HIFUSCORE_DATAMODELS_SONICATION_HPP_INCLUDE

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
  const float mXPos;

  /**
   * @brief The y coordinate for the sonication
   * 
   */
  const float mYPos;

  /**
   * @brief For how long the sonication runs
   * 
   */
  const float mTimeOn;

  /**
   * @brief For how long to cool down after the sonication
   * 
   */
  const float mTimeOff;

  /**
   * @brief The step size for the simulation
   * 
   */
  const float mDT;

  /**
   * @brief Constructs a new auxilliary Sonication object
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
  unsigned getOnSteps() const;

  /**
   * @brief Get the discrete amount of steps for how long the medium cools down after sonication
   * 
   * @return int    - the number of steps
   */
  unsigned getOffSteps() const;

  /**
   * @brief Get the discrete x coordinate for the sonication in the matrix
   * 
   * @return int    - the x coordinate
   */
  unsigned getXPos() const;

  /**
   * @brief Get the discrete y coordinate for the sonication in the matrix
   * 
   * @return int    - the y coordinate
   */
  unsigned getYPos() const;

};
//END of Sonication struct


#endif