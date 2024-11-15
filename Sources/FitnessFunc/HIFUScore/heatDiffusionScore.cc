/**
 * @file        heatDiffusionScore.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The scoreing function for heat diffusion simulation fitness.
 * 
 * @version     0.1
 * 
 * @date        2020-04-06 (created) \n
 */

#include <vector>
#include <cmath>

#include "discreteMap.h"
#include "utils.h"
#include "getHeatSource.h"
#include "calculateHeating.h"
#include "mapsFactory.h"

#include "../FitnessFunction.h"

#ifndef SEP
  #ifdef _WIN32
    #define SEP '\\';
  #else 
    #define SEP '/'
  #endif
#endif

double fitnessFunction(const unsigned n, 
                       double* chromosome)
{
  //std::cout << "fitness called\n" << std::flush;
  const int sonicationsAmount = n / 4;
  const float dt = 0.1f;

  float xPos, yPos, timeOn, timeOff;

  std::vector<Sonication> sonications;

  // to support exploration of optimizers using floating point values we shrink the search space
  for (int i = 0; i < n; i+=4)
  {
    xPos      = chromosome[i]  ;
    yPos      = chromosome[i+1];
    timeOn    = chromosome[i+2];
    timeOff   = chromosome[i+3];

    sonications.push_back(
      Sonication{  
        std::round(xPos),  
        std::round(yPos),  
        std::round(timeOn),  
        std::round(timeOff), 
        dt 
      }
    );
  }

  const std::string dataFolderPath          = std::string("FitnessFunc") + SEP + "HIFUScore" + SEP + "Data";
  const std::string austinWomanFolderPath   = dataFolderPath + SEP + "AustinWoman";
  const std::string targetMapsFolderPath    = dataFolderPath + SEP + "Targets";

  const std::string targetMapCirclesFilePath      = targetMapsFolderPath + SEP + "targetMapCircles.dat";
  const std::string prohibitedMapCirclesFilePath  = targetMapsFolderPath + SEP + "prohibitedMapCircles.dat";

  // static variables to save computation time on multiple fitness executions
  static Medium mediumBeforeTreatment = Medium(austinWomanFolderPath);

  const unsigned nx = mediumBeforeTreatment.nx;
  const unsigned ny = mediumBeforeTreatment.ny;
  const float    dx = mediumBeforeTreatment.dx;

  static DiscreteMap targetMap = GetCirclesMap(nx, 
                                               ny, 
                                               dx, 
                                               targetMapCirclesFilePath,
                                               [](int a, int b) -> int { return a > b ? a + b : b;});

  static DiscreteMap penalizeMap = GetCirclesMap(nx, 
                                                 ny, 
                                                 dx, 
                                                 prohibitedMapCirclesFilePath,
                                                 [](int a, int b) -> int {return a + b;});
  //targetMap.PrintMap(";");
  //std::cout << "\n\n\n";
  //penalizeMap.PrintMap(";");

  Medium medium(mediumBeforeTreatment);
  std::vector<std::vector<float>> heatSources(sonicationsAmount, std::vector<float>(nx*ny, 0.0f));

  #pragma omp parallel for
  for (int i = 0; i < sonicationsAmount; i++)
  {
    heatSources[i] = getHeatSource(medium.nx, medium.ny, sonications[i].xPos, sonications[i].yPos, medium.dx, 1);
  }

  int stillAlive = 0;
  int mistreated = 0;
  for (int i = 0; i < sonicationsAmount; i++)
  {

    medium.q = heatSources[i];
    calculateHeating(medium, sonications[i]);
  }

  DiscreteMap lesionMap = DiscreteMap::CreateLesionMap(medium);
  stillAlive        = (targetMap - (targetMap * lesionMap)).sum();
  mistreated        = (penalizeMap * lesionMap).sum();

  //std::cout << "\n\n\n";
  //lesionMap.PrintMap(";");
  //std::cerr << stillAlive + mistreated << "\n";
  return stillAlive + mistreated;
}


bool isConverged(const unsigned genNumber, 
                 const unsigned fitnessEvaluations,
                 const double optimum,
                 const double bestFitness,
                 const unsigned populationSize,
                 double** populationFitness)
{
  return std::abs(bestFitness - optimum) < 10.0 || fitnessEvaluations > 2000;
}


int _isInConstraints_blob(const unsigned chromoLen, 
                          double *chromozome)
{

  double *upperBounds = new double[4];
  double *lowerBounds = new double[4];
  getConstraints(4, &upperBounds, &lowerBounds);

  double x, y, t_on, t_off;
  int retValue = -1;
  for (int i = 0; i < chromoLen; i+=4)
  {
    x     = chromozome[i];
    y     = chromozome[i+1];
    t_on  = chromozome[i+2];
    t_off = chromozome[i+3];

    if (x > upperBounds[0] || x < lowerBounds[0])
    {
      retValue = i;
      break;
    }
    if (y > upperBounds[1] || y < lowerBounds[1])
    {
      retValue = i + 1;
      break;
    }
    if (t_on > upperBounds[2] || t_on < lowerBounds[2])
    {
      retValue = i + 2;
      break;
    }
    if (t_off > upperBounds[3] || t_off < lowerBounds[3])
    {
      retValue = i + 3;
      break;
    }
  }
  delete[] upperBounds;
  delete[] lowerBounds;
  return retValue;
}

int _isInConstraints_flower(const unsigned chromoLen,
                            double *chromozome)
{
  double *upperBounds = new double[4];
  double *lowerBounds = new double[4];
  getConstraints(4, &upperBounds, &lowerBounds);

  const double innerCircleRadius    = 10;
  const double innerCircleCenter[]  = {247.5, 247.5};
  const double outerCircleRadius    = 40;

  double pointDistanceSquared;
  double radiusSquared;

  //std::cerr << "isInConstraints ";
  double x, y, t_on, t_off;
  int retValue = -1;
  for (int i = 0; i < chromoLen; i+=4)
  {
    x     = chromozome[i];
    y     = chromozome[i+1];
    t_on  = chromozome[i+2];
    t_off = chromozome[i+3];

	//cirlce with cirlce inside
    /*
    double pointDistanceSquared = ((x - innerCircleCenter[0]) * (x - innerCircleCenter[0]) + 
                                   (y - innerCircleCenter[1]) * (y - innerCircleCenter[1]));
    double radiusSquared        = (outerCircleRadius * outerCircleRadius);
    if (pointDistanceSquared > radiusSquared)
    {
      const double distX = std::abs(innerCircleCenter[0] - x);
      const double distY = std::abs(innerCircleCenter[1] - y);

      if (distX < distY)
      {
        //std::cerr << "y" << "\n";
        retValue = i + 1;
        break;
      }
      else
      {
        //std::cerr << "x" << "\n";
        retValue = i;
        break;
      }
    }*/

    if (x > upperBounds[0] || x < lowerBounds[0])
    {
      //std::cerr << "x" << "\n";
      retValue = i;
      break;
    }

    if (y > upperBounds[1] || y < lowerBounds[1])
    {
      //std::cerr << "y" << "\n";
      retValue = i + 1;
      break;
    }

    if (t_on > upperBounds[2] || t_on < lowerBounds[2])
    {
      //std::cerr << "tOn" << "\n";
      retValue = i + 2;
      break;
    }
    if (t_off > upperBounds[3] || t_off < lowerBounds[3])
    {
      //std::cerr << "tOff" << "\n";
      retValue = i + 3;
      break;
    }
    
	// is inside the middle circle, if so, return the axis closest to the middle
    pointDistanceSquared = ((x - innerCircleCenter[0]) * (x - innerCircleCenter[0]) + 
                                   (y - innerCircleCenter[1]) * (y - innerCircleCenter[1]));
    radiusSquared        = (innerCircleRadius * innerCircleRadius);
    if (pointDistanceSquared < radiusSquared)
    {
      const double distX = std::abs(innerCircleCenter[0] - x);
      const double distY = std::abs(innerCircleCenter[1] - y);

      if (distX < distY)
      {
        retValue = i;

      }
      else
      {
        retValue = i + 1;
      }
      break;
    }
    
  }
  delete[] lowerBounds;
  delete[] upperBounds;
  return retValue;
}

/**
 * @brief User defined function to checks if the given chromozome is inside the constraints.
 * 
 * @param chromoLen     - length of the chromozome to check
 * @param chromozome    - the chromozome data 
 * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
 */
int isInConstraints(const unsigned chromoLen, 
                   double* chromozome)
{
  return _isInConstraints_flower(chromoLen, chromozome);
}
// END of isInConstrains

void _applyConstraints_blob(const unsigned n,
                            double** data)
{
  double *upperBounds = new double[4];
  double *lowerBounds = new double[4];
  getConstraints(4, &upperBounds, &lowerBounds);

  double ranges[] = { upperBounds[0] - lowerBounds[0], 
                      upperBounds[1] - lowerBounds[1], 
                      upperBounds[2] - lowerBounds[2], 
                      upperBounds[3] - lowerBounds[3] };

  double x, y, t_on, t_off;
  for (int i = 0; i < n; i += 4)
  {
    x     = data[0][i];
    y     = data[0][i+1];
    t_on  = data[0][i+2];
    t_off = data[0][i+3];

    while (x > upperBounds[0] || x < lowerBounds[0])
    {
      if (x > upperBounds[0])
      { 
        data[0][i] = upperBounds[0] - ((int)std::round(std::abs(x - upperBounds[0])) % (int)ranges[0]);
      }
      else if (x < lowerBounds[0])
      {
        data[0][i] = lowerBounds[0] + ((int)std::round(std::abs(lowerBounds[0] - x)) % (int)ranges[0]);
      }
      x = data[0][i];
    }
    
    while (y > upperBounds[1] || y < lowerBounds[1])
    {
      if (y > upperBounds[1])
      {
        data[0][i+1] = upperBounds[1] - ((int)std::round(std::abs(y - upperBounds[1])) % (int)ranges[1]);
      }
      else if (y < lowerBounds[1])
      {
        data[0][i+1] = lowerBounds[1] + ((int)std::round(std::abs(lowerBounds[1] - y)) % (int)ranges[1]);
      }
      y = data[0][i+1];
    }

    while (t_on > upperBounds[2] || t_on < lowerBounds[2])
    {
      if (t_on > upperBounds[2])
      {
        data[0][i+2] = upperBounds[2] - ((int)std::round(std::abs(t_on - upperBounds[2])) % (int)ranges[2]);
      }
      else if (t_on < lowerBounds[2])
      {
        data[0][i+2] = lowerBounds[2] + ((int)std::round(std::abs(t_on - lowerBounds[2])) % (int)ranges[2]);
      }
      t_on  = data[0][i+2];
    }

    while (t_off > upperBounds[3] || t_off < lowerBounds[3])
    {
      if (t_off > upperBounds[3])
      {
        data[0][i+3] = upperBounds[3] - ((int)std::round(std::abs(t_off - upperBounds[3])) % (int)ranges[3]);
      }
      else if (t_off < lowerBounds[3])
      {
        data[0][i+3] = lowerBounds[3] + ((int)std::round(std::abs(t_off - lowerBounds[3])) % (int)ranges[3]);
      }
      t_off = data[0][i+3];
    }
  }
  delete[] upperBounds;
  delete[] lowerBounds;
}

void _applyConstraints_flower(const unsigned n, 
                              double** data)
{

  double *upperBounds = new double[4];
  double *lowerBounds = new double[4];
  getConstraints(4, &upperBounds, &lowerBounds);

  //const double upperBounds[] = {287.5, 287.5,  5.0, 20.0};
  //const double lowerBounds[] = {207.5, 207.5,  0.0,  1.0};
  double ranges[] = { upperBounds[0] - lowerBounds[0], 
                      upperBounds[1] - lowerBounds[1], 
                      upperBounds[2] - lowerBounds[2], 
                      upperBounds[3] - lowerBounds[3] };

  const double innerCircleRadius    = 10;
  const double innerCircleCenter[]  = {247.5, 247.5};
  
  const double outerCircleRadius    = 40;
  //const double outerCircleCenter[]  = {247.5, 247.5};

  double x, y, t_on, t_off;
  for (int i = 0; i < n; i += 4)
  {
    x     = data[0][i];
    y     = data[0][i+1];
    t_on  = data[0][i+2];
    t_off = data[0][i+3];
  
  /*
    // pushes back inside the circle in the direction created by the line joining the point and the middle 
	// using the distance as a scale.
    double pointDistanceSquared = ((x - innerCircleCenter[0]) * (x - innerCircleCenter[0]) + 
                                   (y - innerCircleCenter[1]) * (y - innerCircleCenter[1]));
    double radiusSquared        = (outerCircleRadius * outerCircleRadius);
    if (pointDistanceSquared > radiusSquared)
    {
      //std::cerr << " = ";
      double scale = (std::sqrt(pointDistanceSquared)/outerCircleRadius) - 1.0;
      double shift = outerCircleRadius * scale;

      double theta = std::atan2(y - innerCircleCenter[1], x - innerCircleCenter[0]);
      double intersection[] = { 
          innerCircleCenter[0] + (outerCircleRadius - shift) * std::cos(theta), 
          innerCircleCenter[1] + (outerCircleRadius - shift) * std::sin(theta)
      };
      x = data[0][i]    = intersection[0];
      y = data[0][i+1]  = intersection[1];
      //std::cerr << x << " " << y << "\n";
    }
  */

    while (x > upperBounds[0] || x < lowerBounds[0])
    {
      //std::cerr << "tOn apply constraints " << t_on << "\n";
      if (x > upperBounds[0])
      {
        data[0][i] = upperBounds[0] - ((int)std::round(std::abs(x - upperBounds[0])) % (int)ranges[0]);
      }
      else if (x < lowerBounds[0])
      {
        data[0][i] = lowerBounds[0] + ((int)std::round(std::abs(x - lowerBounds[0])) % (int)ranges[0]);
      }
      x  = data[0][i];
    }

    while (y > upperBounds[1] || y < lowerBounds[1])
    {
      //std::cerr << "tOff apply constraints " << t_off << "\n";
      if (y > upperBounds[1])
      {
        data[0][i+1] = upperBounds[1] - ((int)std::round(std::abs(y - upperBounds[1])) % (int)ranges[1]);
      }
      else if (y < lowerBounds[1])
      {
        data[0][i+1] = lowerBounds[1] + ((int)std::round(std::abs(y - lowerBounds[1])) % (int)ranges[1]);
      }
      y = data[0][i+1];
    }

    while (t_on > upperBounds[2] || t_on < lowerBounds[2])
    {
      //std::cerr << "tOn apply constraints " << t_on << "\n";
      if (t_on > upperBounds[2])
      {
        data[0][i+2] = upperBounds[2] - ((int)std::round(std::abs(t_on - upperBounds[2])) % (int)ranges[2]);
      }
      else if (t_on < lowerBounds[2])
      {
        data[0][i+2] = lowerBounds[2] + ((int)std::round(std::abs(t_on - lowerBounds[2])) % (int)ranges[2]);
      }
      t_on  = data[0][i+2];
    }

    while (t_off > upperBounds[3] || t_off < lowerBounds[3])
    {
      //std::cerr << "tOff apply constraints " << t_off << "\n";
      if (t_off > upperBounds[3])
      {
        data[0][i+3] = upperBounds[3] - ((int)std::round(std::abs(t_off - upperBounds[3])) % (int)ranges[3]);
      }
      else if (t_off < lowerBounds[3])
      {
        data[0][i+3] = lowerBounds[3] + ((int)std::round(std::abs(t_off - lowerBounds[3])) % (int)ranges[3]);
      }
      t_off = data[0][i+3];
    }
  
  
  /*
    x     = data[0][i];
    y     = data[0][i+1];

    double pointDistanceSquared = ((x - innerCircleCenter[0]) * (x - innerCircleCenter[0]) + (y - innerCircleCenter[1]) * (y - innerCircleCenter[1]));
    double radiusSquared        = (innerCircleRadius * innerCircleRadius);
    if (pointDistanceSquared < radiusSquared)
    {
      //std::cerr << "applyConstraints " << x << " " << y << " = ";
      double scale = std::sqrt(pointDistanceSquared)/innerCircleRadius;
      double shift = (outerCircleRadius - innerCircleRadius) * scale;

      double theta = std::atan2(y - innerCircleCenter[1], x - innerCircleCenter[0]);
      double intersection[] = { 
          innerCircleCenter[0] + (innerCircleRadius + shift) * std::cos(theta), 
          innerCircleCenter[1] + (innerCircleRadius + shift) * std::sin(theta)
      };
      x = data[0][i]    = intersection[0];
      y = data[0][i+1]  = intersection[1];
      //std::cerr << "" << x << " " << y << "\n";
    }
  */
  
  }
  delete[] upperBounds;
  delete[] lowerBounds;
}

/**
 * @brief user defined way of how to handle the current optimized variables that are outside of constraints
 * 
 * @param [in] n          - size of the data (for problem dimension 'n' and population size 'p' its 'n * p')
 * @param [in, out] data  - current data of points in search space (variables) 
 */
void applyConstraints(const unsigned n, 
                      double** data)
{
  _applyConstraints_flower(n, data);
}
// END of applyConstraints

void _getConstraints_flower(const unsigned n, 
                            double** uB,
                            double** lB)
{
  const double upperBounds[] = {287.5, 287.5, 10.0, 120.0};
  const double lowerBounds[] = {207.5, 207.5,  1.0,   1.0};

  for (int i = 0; i < n; i += 4)
  {
    uB[0][i] = upperBounds[0];
    lB[0][i] = lowerBounds[0];

    uB[0][i+1] = upperBounds[1];
    lB[0][i+1] = lowerBounds[1];

    uB[0][i+2] = upperBounds[2];
    lB[0][i+2] = lowerBounds[2];

    uB[0][i+3] = upperBounds[3];
    lB[0][i+3] = lowerBounds[3];
  }
}

void _getConstraints_blob(const unsigned n, 
                          double** uB,
                          double** lB)
{
  const double upperBounds[] = {342.0, 293.0, 20.0, 60.0};
  const double lowerBounds[] = {272.0, 233.0,  1.0,  1.0};
  for (int i = 0; i < n; i += 4)
  {
    uB[0][i] = upperBounds[0];
    lB[0][i] = lowerBounds[0];

    uB[0][i+1] = upperBounds[1];
    lB[0][i+1] = lowerBounds[1];

    uB[0][i+2] = upperBounds[2];
    lB[0][i+2] = lowerBounds[2];

    uB[0][i+3] = upperBounds[3];
    lB[0][i+3] = lowerBounds[3];
  }
}

/**
 * @brief way for the user to set outer bounds on searched variables (domain constraints)
 * 
 * @param [in] n                  - dimension of the problem (number of searched variables)
 * @param [in, out] upperBounds   - upper constraints on search space for each variable
 * @param [in, out] lowerBounds   - lower constraints on search space for each variable
 */
void getConstraints(const unsigned n, 
                    double** uB,
                    double** lB)
{
  _getConstraints_flower(n, uB, lB);
}
//END of getConstraints 