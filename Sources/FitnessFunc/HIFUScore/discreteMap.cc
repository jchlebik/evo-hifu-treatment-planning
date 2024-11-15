/**
 * @file        discreteMap.c
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       A source code file containing all definitions for DiscreteMap class.
 * 
 * @version     0.2
 * 
 * @date        2020-03-12 (created) \n
 *              2020-03-18 (update) \n
 */

#include "discreteMap.h"

/**
 * @brief Construct a new 2D TargetMap object of xDimSize to yDimSize with circles specified in coordinates. \n
  *        Edge cases are wrapped around its opposite edge calculated by modulus operation. Identical to the 
  *        Matlab version used for fitness maps in k-Wave EUD.
  * 
  * @param [in] g              - a Grid class carrying relevant information about the map (xn, yn, dx)
  * @param [in] data           - the flattened data values of the map.
  */  
DiscreteMap::DiscreteMap(GridInfo g,
                         std::vector<int>& mapData):
GridInfo(g), data(mapData)
{}
//END of TargetMap(Grid) ctor

/**
 * @brief Construct a new 2D TargetMap object of xDimSize to yDimSize with circles specified in coordinates. \n
  *       Edge cases are wrapped around its opposite edge calculated by modulus operation. Result is  \n 
  *       dentical to the Matlab version used for fitness maps in k-Wave HIFU optimization.
  * 
  * @param [in] xDimSize       - the size of the x axis of the discrete medium space
  * @param [in] yDimSize       - the size of the y axis of the discrete medium space
  * @param [in] dx             - the distance between two points in the discrete space in meters
  * @param [in] data           - the flattened data values of the map.
  */
DiscreteMap::DiscreteMap(const unsigned xDimSize, 
                         const unsigned yDimSize, 
                         const float dx,
                         std::vector<int>& mapData): 
GridInfo(int(xDimSize), int(yDimSize), dx), data(mapData)
{}
// END of DiscreteMap(unsigned, unsigned, float) ctor

/**
 * @brief sums the data contained in the map
 * 
 * @return int - sum of the data in the map
 */
int DiscreteMap::sum()
{
  int sumOfElements = 0;

  //#pragma omp parallel for reduction(+:sumOfElements)
  for (int i = 0; i < data.size(); i++)
  {
    sumOfElements += data[i];
  }
  //std::for_each(data.begin(), data.end(), [&sumOfElements] (int n) { sumOfElements += n; });
  return sumOfElements;
}
//END of sum

/**
 * @brief Definition of the multiplication operation of the discrete maps. \n
 *        Equivalent to scalar product of the matrices. Requires L an R to \n
 *        be of the same size
 * 
 * @param [in] L          - L side argument of the multiplication
 * @param [in] R          - R side argument of the multiplication
 * @return DiscreteMap    - L * R
 */
DiscreteMap operator*(DiscreteMap const &L, DiscreteMap const &R)
{
  if (L.data.size() != R.data.size()) 
  { 
    throw std::logic_error("Incompatible map sizes for operator."); 
  }

  std::vector<int> res(L.data.size());
  for (int i = 0; i < res.size(); i++)
  {
    res[i] = L.data[i] * R.data[i];
  }
  //std::transform(L.data.begin(), L.data.end(), R.data.begin(), res.begin(), [](const int &a, const int &b) {return a * b;});
  return DiscreteMap(L.nx, L.ny, L.dx, res);
}
//END of operator*()

/**
 * @brief Definition of the substration operation for the discrete maps. \n
 *        Equivalent to substraction of the matrices. Requires L an R to \n
 *        be of the same size
 * 
 * @param [in] L          - L side argument of the substraction
 * @param [in] R          - R side argument of the substraction
 * @return DiscreteMap    - L - R
 */
DiscreteMap operator-(DiscreteMap const &L, DiscreteMap const &R)
{
  if (L.data.size() != R.data.size()) 
  { 
    throw std::logic_error("Incompatible map sizes for operator.");
  }

  std::vector<int> res(L.data.size());
  for (int i = 0; i < res.size(); i++)
  {
    res[i] = L.data[i] - R.data[i];
  }
  //std::transform(L.data.begin(), L.data.end(), R.data.begin(), res.begin(), [](const int &a, const int &b) {return a - b;});
  return DiscreteMap(L.nx, L.ny, L.dx, res);
}
//END of operator-()

/**
 * @brief Definition of the addition operation of the discrete maps. \n
 *        Equivalent to dot addition of matrices. Requires L an R to \n
 *        be of the same size
 * 
 * @param [in] L          - L side argument of the addition
 * @param [in] R          - R side argument of the addition
 * @return DiscreteMap    - L .+ R
 */
DiscreteMap operator+(DiscreteMap const &L, DiscreteMap const &R)
{
  if (L.data.size() != R.data.size()) 
  { 
    throw std::logic_error("Incompatible map sizes for operator.");
  }

  std::vector<int> res(L.data.size());
  for (int i = 0; i < res.size(); i++)
  {
    res[i] = L.data[i] + R.data[i];
  }
  //std::transform(L.data.begin(), L.data.end(), R.data.begin(), res.begin(), [](const int &a, const int &b) {return a - b;});
  return DiscreteMap(L.nx, L.ny, L.dx, res);
}
//END of operator+()


/**
 * @brief Creates a lesion map specific to the current state of the given medium
 * 
 * @param [in] m             - Medium instance to create the lesion map for
 * @return DiscreteMap       - a binary DiscreteMap representing the lesion map
 */
DiscreteMap DiscreteMap::CreateLesionMap(Medium& m)
{
  std::vector<int> lesionMap;
  lesionMap.reserve(m.cem43.size());

  for (int i = 0; i < m.cem43.size(); i++)
  {
    if (m.cem43[i] >= 240)
    {
      lesionMap.push_back(1);
    }
    else
    {
      lesionMap.push_back(0);
    }
  }
  //auto filter = [&lesionMap](int val) -> void { (val > 240) ? lesionMap.push_back(1) : lesionMap.push_back(0); };
  //std::for_each(m.cem43.begin(), m.cem43.end(), filter);
  return DiscreteMap(m.nx, m.ny, m.dx, lesionMap);
}
// END of CreateLesionMap

/**
 * @brief Creates a DiscreteMap object of the required size from the circles defining file
 * 
 * @param [in] nx                     - size of the map on x dimension
 * @param [in] ny                     - size of the map on y dimension
 * @param [in] dx                     - distance between the neighbours in the real medium [m]
 * @param [in] circles                - a vector container of the CircleInfo instances representing the circles on the \n
 *                                      discrete map
 * @param [in] ftor                   - a function - how apply the value to its position (e.g. by addition or replacing)
 * @return DiscreteMap                - a nx x ny map filled with circles given in container
 */
DiscreteMap DiscreteMap::CreateTargetsMap(const unsigned nx, 
                                          const unsigned ny, 
                                          const float dx, 
                                          std::vector<CircleInfo>& circles,
                                          std::function<int(int, int)> ftor)
{

  //int *map = (int*)aligned_alloc(64, nx * ny * sizeof(int));
  std::vector<int> res(nx*ny, 0);

  //std::vector<int> map(nx * ny, 0);

  // negative modulo lambda function
  auto mod = [](int a, int b) -> int{
    if (b == 0)
    {
      return (unsigned)!((int)0);
    }
    int ret = a % b;
    while(ret < 0)
    {
      ret += b;
    }
    return ret; 
  };
  // instead of iterating over the entire matrix and calculating distances to every 
  // given center, we create a "box" around each circle and 
  // iterate over fields contained in that box
  for(auto &coord : circles)
  {
    // find box bounds
    // with edge cases, the bounds are allowed to be negative for smooth iterations (extending the marix)
    int lowBound    = coord.x - coord.radius;
    int highBound   = coord.x + coord.radius;
    int leftBound   = coord.y - coord.radius;
    int rightBound  = coord.y + coord.radius;

    // iterate over the fields in box
    //#pragma omp parallel for
    for (int row = lowBound; row < highBound; row++)
    {
      int moduloRow = mod(row, int(nx)); // cmath modulus does not calculate negative modulos correctly
      //#pragma omp simd
      for (int col = leftBound; col < rightBound; col++)
      {
        int moduloCol = mod(col, int(ny));
        int index = moduloRow * nx + moduloCol;

        // no need to round the distance because we are using <. we can just cut the floating part off using integer
        //int distance = std::sqrt(std::pow((row - (int(coord.x) - 1)), 2) + std::pow((col - (int(coord.y) - 1)), 2));
        //if (distance < coord.radius)
        //{
        //  map[index] += coord.circleValue; // overlapping values are added
        //}

        // faster than calculating the square root of distance and power operation is congruent on Inequality split
        int distanceSquared = (row - (int(coord.x) - 1)) * (row - (int(coord.x) - 1))  + (col - (int(coord.y) - 1)) * (col - (int(coord.y) - 1));
        int radiusSquared   = coord.radius * coord.radius;  

        if (distanceSquared < radiusSquared)
        {
          res[index] = ftor(res[index], coord.circleValue);
        }
      }
    }
    //std::vector<int> res(map, map + nx*ny);
    //free(map);
  }
  return DiscreteMap(nx, ny, dx, res);
}
  //END of CreateTargetsMap
