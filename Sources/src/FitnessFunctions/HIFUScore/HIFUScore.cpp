#include <FitnessFunctions/FitnessFunction.hpp>
#include <FitnessFunctions/HIFUScore.hpp>

#include <FitnessFunctions/HIFUScore/DataModels/GridInfo.hpp>
#include <FitnessFunctions/HIFUScore/DataModels/DiscreteMap.hpp>
#include <FitnessFunctions/HIFUScore/DataModels/Medium.hpp>
#include <FitnessFunctions/HIFUScore/DataModels/Sonication.hpp>
#include <FitnessFunctions/HIFUScore/DataModels/utils.hpp>

#include <FitnessFunctions/HIFUScore/KWaveDiffusionSolver/kWaveDiffusionSolver.hpp>

#include <sstream>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <limits>
#include <vector>
#include <string>

#include <cmath>

/**
 * @brief Create a Gaussian distribution on specified base with specified parameters.
 * 
 * @param [in] startBase          - the starting point on the axis
 * @param [in] endBase            - the ending point on the axis
 * @param [in] magnitude          - bell height (defaul is normalised)
 * @param [in] mean               - mean or expected value (default = 0)
 * @param [in] variance           - bell width (default = 1)
 * @return std::vector<float>     - a Gaussian distribution $f(base)$ with the specified magnitude, mean and variance
 */
std::vector<float> gaussianDistribution(const float startBase,
                                        const float endBase,
                                        float magnitude = std::numeric_limits<float>::max(),
                                        float mean = 0.0,
                                        float variance = 1.0)
{

  std::vector<float> gauss;
  if (magnitude == std::numeric_limits<float>::max())
  {
    magnitude = std::pow((2 * M_PI * variance), -0.5);    // normalize in case magnitude is not given
  }

  // save some calculations by taking out invariant division
  const float divCoeff = 1.0f / (2.0f * variance);
  float currentBaseValue = startBase;

  while (currentBaseValue <= endBase)
  {
    float xMinusMean = currentBaseValue - mean;  // x * x shows faster results than pow(x, 2)
    float value = magnitude * std::exp( -(xMinusMean * xMinusMean) * divCoeff ); 
    gauss.push_back(value);
    currentBaseValue += 1.0f;
  }

  return gauss;
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
float* getHeatSource(const unsigned nx,
                      const unsigned ny,
                      const unsigned x,
                      const unsigned y,
                      const float    dx,
                      const float    magnitude)
{
  //constants taken from the matlab script of the same name
  const float weightCoeff       = 10e6;
  //const float magnitude         = 1;
  const float fullWidthHalfMax  = 10e-3;
  const float variance          = std::pow((fullWidthHalfMax / (2.0f * dx)), 2) / (2.0f * std::log(2.0f));

  std::vector<float> gaussX = gaussianDistribution(1.0f, static_cast<float>(nx), magnitude, x, variance / 8.0f);
  std::vector<float> gaussY = gaussianDistribution(1.0f, static_cast<float>(ny), magnitude, y, variance);

  const unsigned rowsNum = gaussY.size();
  const unsigned colsNum = gaussX.size();

  float* heatMatrix = reinterpret_cast<float*>(_mm_malloc(sizeof(float) * rowsNum * colsNum, _ALIGN_)); //flattened

  // here we need to do gaussX^T * gaussY
  // because we are dealing with (1 x Nx) and (1 x Ny) matrixes and 
  // (Nx x 1)*(1 x Ny) matrix multiplication, we can simplify the operation to
  // O(n^2)
  #pragma omp parallel for schedule(guided)
  for (unsigned row = 0; row < rowsNum; row++)
  {
    float value             = gaussX[row];
    unsigned oneDIndexStart = row * rowsNum;
    
    #pragma omp simd
    for (unsigned col = 0; col < colsNum; col++)
    {
      heatMatrix[oneDIndexStart + col] = weightCoeff * value * gaussY[col];
    }
  }
  return heatMatrix;
}
// end of getHeatSource



HIFUScore::HIFUScore(ProblemType problemType, 
                     const unsigned numOfSonications,
                     const float optValue,
                     const float epsilon):
  FitnessFunction(optValue, epsilon),
  mProblemType(problemType),
  mNumberOfSonications(numOfSonications)
{
  const std::string dataFolderPath          = "Data";
  const std::string austinWomanFolderPath   = dataFolderPath + SEP + "AustinWoman";
  const std::string targetMapsFolderPath    = dataFolderPath + SEP + "Targets";

  std::string problemTypeString;
  if(problemType == ProblemType::Flower)
  {
    problemTypeString = std::string("Flower");
  }
  else
  {
    problemTypeString = std::string("Blob");
  }

  const std::string targetMapCirclesFilePath      = targetMapsFolderPath + SEP + problemTypeString + SEP + "targetMapCircles.dat";
  const std::string prohibitedMapCirclesFilePath  = targetMapsFolderPath + SEP + problemTypeString + SEP + "prohibitedMapCircles.dat";
  const std::string constraintsFilePath           = targetMapsFolderPath + SEP + problemTypeString + SEP + "constraints.dat";

  mMediumBeforeTreatment.reset(new Medium(austinWomanFolderPath));

  mTargetMap.reset(
    DiscreteMap::getCirclesMap(mMediumBeforeTreatment->mNx,
                               mMediumBeforeTreatment->mNy,
                               mMediumBeforeTreatment->mDx,
                               targetMapCirclesFilePath,
                               [](int a, int b) -> int { return a > b ? a + b : b;})
  );

  mPenalizeMap.reset(
    DiscreteMap::getCirclesMap(mMediumBeforeTreatment->mNx,
                               mMediumBeforeTreatment->mNy,
                               mMediumBeforeTreatment->mDx,
                               prohibitedMapCirclesFilePath,
                               [](int a, int b) -> int {return a + b;})
  );

  mKWaveDiffSolver.reset(
    new kWaveDiffusionSolver2D(mMediumBeforeTreatment->mNx,     // medium grid size on x axis
                               mMediumBeforeTreatment->mNy,     // medium grid size on y axis
                               mMediumBeforeTreatment->mRho,    // medium density matrix - tissue mass density [kg/m^3]
                               mMediumBeforeTreatment->mC,      // medium specific heat capacity matrix - tissue specific heat capacity [J/(kg.K)]
                               mMediumBeforeTreatment->mLambda) // medium thermal conductivity - tissue thermal conductivity [W/(m.K)]
  );

  setupConstraints(constraintsFilePath);
}

float HIFUScore::operator()(const std::vector<float>& chromosome)
{
  const float dt = 0.1f;
  float xPos, yPos, timeOn, timeOff;

  std::vector<Sonication> sonications;

  for (int i = 0; i < chromosome.size(); i += 4)
  {
    xPos      = chromosome[i]  ;
    yPos      = chromosome[i+1];
    timeOn    = chromosome[i+2];
    timeOff   = chromosome[i+3];

    sonications.push_back(
      { std::round(xPos), std::round(yPos), std::round(timeOn), std::round(timeOff), dt }
    );
  }
  
  const int magnitude = 1;
  std::vector<float*> heatSources(sonications.size());

  for (int i = 0; i < sonications.size(); i++)
  {
    heatSources[i] = getHeatSource(mMediumBeforeTreatment->mNx, 
                                   mMediumBeforeTreatment->mNy, 
                                   sonications[i].mXPos, 
                                   sonications[i].mYPos, 
                                   mMediumBeforeTreatment->mDx, 
                                   magnitude);
  }

  Medium mediumCopy(*mMediumBeforeTreatment);
  for (int i = 0; i < sonications.size(); i++)
  {
    mediumCopy.setHeatDeposit(heatSources[i]);
    mKWaveDiffSolver->setup(mediumCopy.mQ, mediumCopy.mT, mediumCopy.mCem43);

    // perform the sonication
    int onSteps = sonications[i].getOnSteps();
    if (onSteps != 0)
    {
      mKWaveDiffSolver->takeTimeStep(onSteps,  true,  sonications[i].mDT); // perform heating
    }

    // cooling interval
    std::fill(mediumCopy.mQ, mediumCopy.mQ + mediumCopy.mMatSize, 0.0f);
    int offSteps = sonications[i].getOffSteps();
    if (offSteps != 0)
    {
      mKWaveDiffSolver->takeTimeStep(offSteps, false, sonications[i].mDT); // cool down
    }

    _mm_free(heatSources[i]);
  }

  DiscreteMap lesionMap = DiscreteMap::createLesionMap(mediumCopy);
  int stillAlive        = (*mTargetMap - (*mTargetMap * lesionMap)).sum();
  int mistreated        = (*mPenalizeMap * lesionMap).sum();

  return static_cast<float>(stillAlive + mistreated);
}

void HIFUScore::applyConstraints(std::vector<float>& data)
{
  if (mConstraintType == ConstraintType::Box)
  {
    applyConstraintsForBox(data);
  }
  else if (mConstraintType == ConstraintType::Circle)
  {
    applyConstraintsForCircle(data);
  }
}

int HIFUScore::isInConstraints(const std::vector<float>& chromozome)
{
  if (mConstraintType == ConstraintType::Box)
  {
    isInConstraintsForBox(chromozome);
  }
  else if (mConstraintType == ConstraintType::Circle)
  {
    isInConstraintsForCircle(chromozome);
  }
}

void HIFUScore::getConstraints(const unsigned n, 
                               std::vector<float>& upperBounds,    // out var
                               std::vector<float>& lowerBounds)    // out var
{
  if (n % upperBounds.size() != 0)
  {
    throw std::runtime_error("The requested bounds dimension is not divisible by the constraints provided.");
  }

  for (int i = 0; i < n; i += upperBounds.size())
  {
    for (int j = 0; j < upperBounds.size(); j++)
    {
      upperBounds[i + j] = mUpperBounds[j];
      lowerBounds[i + j] = mLowerBounds[j];
    }
  }
}

/**
 * @brief Way of handling the current variables that are outside of bounds. \n
 *        Specialized implementation for the Box constraitns type.
 * 
 * @param [in, out] data  - variables inside the search space to constraint.
 */
void HIFUScore::applyConstraintsForBox(std::vector<float>& data)
{
  for (int i = 0; i < data.size(); i++)
  {
    int boundIndex = i % mUpperBounds.size();
    int range = mUpperBounds[boundIndex] - mLowerBounds[boundIndex];

    int distance = 0;
    while (data[i] > mUpperBounds[boundIndex] || data[i] < mLowerBounds[boundIndex])
    {
      if (data[i] > mUpperBounds[boundIndex])
      { 
        distance = static_cast<int>(std::round(std::abs(data[i] - mUpperBounds[boundIndex])));
        data[i] = mUpperBounds[boundIndex] - (distance % range);
      }
      else if (data[i] < mLowerBounds[boundIndex])
      {
        distance = static_cast<int>(std::round(std::abs(mLowerBounds[boundIndex] - data[i])));
        data[i] = mLowerBounds[boundIndex] + (distance % range);
      }
    }
  }
}

/**
 * @brief Way of handling the current variables that are outside of bounds. \n
 *        Specialized implementation for the Circle constraitns type.
 * 
 * @param [in, out] data  - variables inside the search space to constraint.
 */
void HIFUScore::applyConstraintsForCircle(std::vector<float>& data)
{
  //TODO
  applyConstraintsForBox(data);
}

/**
 * @brief Checks if the given chromozome is inside the constraints without modifying it. \n 
 *        Specialized implementation for the Box constraitns type.
 * 
 * @param chromozome    - the chromozome data 
 * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
 */
int HIFUScore::isInConstraintsForBox(const std::vector<float>& chromozome)
{
  for (int i = 0; i < chromozome.size(); i++)
  {
    int boundIndex = i % mUpperBounds.size();

    if (chromozome[i] > mUpperBounds[boundIndex] || chromozome[i] < mLowerBounds[boundIndex])
    {
      return i;
    }
  }
  return -1;  
}

/**
 * @brief Checks if the given chromozome is inside the constraints without modifying it. \n 
 *        Specialized implementation for the Circle constraitns type.
 * 
 * @param chromozome    - the chromozome data 
 * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
 */
int HIFUScore::isInConstraintsForCircle(const std::vector<float>& chromozome)
{
  return isInConstraintsForBox(chromozome);
}

/**
 * @brief Setups the search space constraints data. All data read from a file. First line in file defines \n 
 *        the type of constraint: 'B' for box and 'C' for circle. Second line are the lower bounds separated \n 
 *        by comma. Third line are the upper bounds separated by comma. For the 'C' type, next line defines the \n
 *        constraining circle in format of 'centerX, centerY, radius'. Bounds are periodically repeated in case \n
 *        where the problem dimension is higher than given number of bounds. Example inside Data/Blob/constraints.dat
 * 
 * @param [in] constraitnsFilePath      - a file defining the constraints. 
 * 
 */
void HIFUScore::setupConstraints(const std::string& constraintsFilePath)
{
  std::stringstream lineStream;
  std::string lineInputString;

  std::fstream constraintsStream(constraintsFilePath, std::ios_base::in);    
  if (constraintsStream.bad())
  {
    throw std::runtime_error("File " + constraintsFilePath + " with constraints specification could not be opened.");
  }

  constraintsStream >> std::ws; 
  std::getline(constraintsStream, lineInputString);
  std::transform(lineInputString.begin(), lineInputString.end(), lineInputString.begin(),
    [](unsigned char c){ 
      return std::tolower(c); 
  });

  if (lineInputString == "s" || lineInputString == "b") //square or box
  {
    mConstraintType = ConstraintType::Box;
  }
  else if (lineInputString == "c") //circular constraint
  {
    mConstraintType = ConstraintType::Circle;
  }
  else
  {
    throw std::runtime_error("Unknown constraints type. Current support are B for box or C for circle.");
  }

  std::vector<std::reference_wrapper<std::vector<float>>> bounds;

  if (mConstraintType == ConstraintType::Box)
  {
    bounds = std::vector<std::reference_wrapper<std::vector<float>>>
    {
      std::ref(mLowerBounds),
      std::ref(mUpperBounds)
    };

    constraintsStream >> std::ws;
    for (int i = 0; i < 2; i++)
    {
      constraintsStream >> std::ws;
      std::getline(constraintsStream, lineInputString);
      lineStream.str(lineInputString);

      std::string val;
      lineStream >> std::ws;
      while (std::getline(lineStream,  val,  ','))
      {
        bounds[i].get().push_back(std::stof(val));
        lineStream >> std::ws;
      }
    }
  }
  else if (mConstraintType == ConstraintType::Circle)
  {
    //TODO
  }
}