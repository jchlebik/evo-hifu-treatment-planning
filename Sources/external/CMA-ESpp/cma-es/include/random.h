#pragma once

#include <ctime>
#include <cmath>

namespace CMAESpp
{
/**
 * @class Random
 * A pseudo random number generator.
 */
template<typename T>
class Random
{
  // variables for uniform()
  long int startseed;
  long int aktseed;
  long int aktrand;
  long int rgrand[32];
  // variables for gauss()
  bool     stored;
  T        hold;
 public:
  /**
   * @param seed use clock if 0
   */
  Random(long unsigned seed = 0);
  /**
   * @param seed 0 == 1
   */
  void start(long unsigned seed);
  /**
   * @return (0,1)-normally distributed random number
   */
  T gauss(void);
  /**
   * @return (0,1)-uniform distributed random number
   */
  T uniform(void);
};
}