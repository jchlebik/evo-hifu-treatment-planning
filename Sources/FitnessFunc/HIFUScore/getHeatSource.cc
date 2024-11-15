/**
 * @file        GetHeatSource.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for creating a Gaussian heating source
 * 
 * @version     0.3
 * 
 * @date        2020-03-02 (created) \n
 *              2020-03-06 (update)
 */

#include "getHeatSource.h"

/**
 * @brief Create a Gaussian distribution on specified base with specified parameters.
 * 
 * @param [in] base               - the x axis on which the distribution is made
 * @param [in] magnitude          - bell height (defaul is normalised)
 * @param [in] mean               - mean or expected value (default = 0)
 * @param [in] variance           - bell width (default = 1)
 * @return std::vector<float>     - a Gaussian distribution $f(base)$ with the specified magnitude, mean and variance
 */
std::vector<float> gaussianDistribution(std::vector<float> base,
                                    float magnitude = std::numeric_limits<float>::max(),
                                    float mean = 0.0,
                                    float variance = 1.0)
{
  if (magnitude == std::numeric_limits<float>::max())
  {
    magnitude = std::pow((2 * M_PI * variance), -0.5);    // normalize in case magnitude is not given
  }

  // save some calculations by taking out invariant division
  const float divCoeff = 1.0f / (2.0f * variance);
  for (float& x : base)
  {
    float xMinusMean = x - mean;  // x * x shows faster results than pow(x, 2)
    float value = magnitude * std::exp( -(xMinusMean * xMinusMean) * divCoeff ); 
    //float value = magnitude * std::exp( -(std::pow((x - mean), 2)) / (2.0f * variance) );
    x = float(value);
  }
  return base;
}
// end of gaussianDistribution

/**
 * @brief Create a Gaussian heating source centered at the position x, y \n
 *        discretised using the given grid parameters.
 * 
 * @param [in] nx     - the amount of points in grid on the x axis
 * @param [in] ny     - the amount of points in grid on the y axis
 * @param [in] x      - x coordinate of the heating source position
 * @param [in] y      - y coordinate of the heating source position
 * @param [in] dx     - spacing of actual samples in the grid (in meters)
 * @return std::vector<float> 
 */
std::vector<float> getHeatSource(const unsigned nx,
                                 const unsigned ny,
                                 const unsigned x,
                                 const unsigned y,
                                 const float dx,
                                 const float magnitude)
{
  //constants taken from the matlab script of the same name
  const float weightCoeff       = 10e6;
  //const float magnitude         = 1;
  const float fullWidthHalfMax  = 10e-3;
  const float variance          = std::pow((fullWidthHalfMax / (2.0f * dx)), 2) / (2.0f * std::log(2.0f));
  
  float n = 0.0;
  auto generator = [n]() mutable {n += 1.0f; return float(n);};

  std::vector<float> xAxis(nx);
  std::vector<float> yAxis(ny);

  // fill the axis and generate gaussian distribution 
  // TODO optimize this code, possible to do in one cycle instead of four
  std::generate(xAxis.begin(), xAxis.end(), generator);
  std::generate(yAxis.begin(), yAxis.end(), generator);
  std::vector<float> gaussX = gaussianDistribution(xAxis, magnitude, x, variance / 8.0f);
  std::vector<float> gaussY = gaussianDistribution(yAxis, magnitude, y, variance);

  std::vector<float> heatMatrix(gaussX.size() * gaussY.size()); //flattened
  const unsigned rowsNum = gaussX.size();
  const unsigned colsNum = gaussY.size();

  // here we need to do gaussX^T * gaussY
  // because we are dealing with (1 x Nx) and (1 x Ny) matrixes and 
  // (Nx x 1)*(1 x Ny) matrix multiplication, we can simplify the operation to
  // O(n^2)
  for (unsigned row = 0; row < rowsNum; row++)
  {
    float value             = gaussX[row];
    unsigned oneDIndexStart = row * rowsNum;
    for (unsigned col = 0; col < colsNum; col++)
    {
        heatMatrix[oneDIndexStart + col] = weightCoeff * value * gaussY[col];
    }
  }
  return heatMatrix;
}
// end of getHeatSource