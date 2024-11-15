#pragma once
#ifndef KEVOOPT_DATAMODELS_MIGRATION_HPP_INCLUDE
#define KEVOOPT_DATAMODELS_MIGRATION_HPP_INCLUDE

#include <vector>

template <class T>
struct Emigrants
{
  std::vector<T> fitness;
  std::vector<T> population;  //flattened

  Emigrants(const size_t numberOfImmigrants, 
            const size_t problemDimension) : 
    fitness(numberOfImmigrants), 
    population(numberOfImmigrants * problemDimension)
  {
  }
};

template <class T>
struct Immigrants
{
  std::vector<T> fitness;
  std::vector<T> population;  //flattened

  Immigrants(const size_t numberOfImmigrants, 
             const size_t problemDimension):
    fitness(numberOfImmigrants),
    population(numberOfImmigrants * problemDimension)
  {}
};
#endif