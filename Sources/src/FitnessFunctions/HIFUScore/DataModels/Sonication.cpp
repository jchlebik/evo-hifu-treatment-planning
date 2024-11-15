/**
 * @file        Sonication.cpp
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Source code file containing all definition for auxilliar classes used during the simulation.
 * 
 * @version     1.0
 * 
 * @date        2020-03-02 (created) \n
 *              2020-02-23 (update) \n
 */

#include <FitnessFunctions/HIFUScore/DataModels/Sonication.hpp>

#include <cmath>

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
  mXPos(xP),
  mYPos(yP),
  mTimeOn(tOn),
  mTimeOff(tOff),
  mDT(dt) {}
// END of Sonication::ctor

/**
 * @brief Get the discrete amount of steps for how long the sonication runs
 * 
 * @return int    - the number of steps
 */
unsigned Sonication::getOnSteps() const
{
  return static_cast<unsigned>(std::round(mTimeOn/mDT));
}
// END of Sonication::getOnSteps

/**
 * @brief Get the discrete amount of steps for how long the medium cools down after sonication
 * 
 * @return int    - the number of steps
 */
unsigned Sonication::getOffSteps() const
{
  return static_cast<unsigned>(std::round(mTimeOff/mDT));
}
// END of Sonication::getOffSteps

/**
 * @brief Get the discrete x coordinate for the sonication in the matrix
 * 
 * @return int    - the x coordinate
 */
unsigned Sonication::getXPos() const
{
  return static_cast<unsigned>(std::round(mXPos));
}
// END of Sonication::getXPos

/**
 * @brief Get the discrete y coordinate for the sonication in the matrix
 * 
 * @return int    - the y coordinate
 */
unsigned Sonication::getYPos() const
{
  return static_cast<unsigned>(std::round(mYPos));
}
// END of Sonication::getYPos
