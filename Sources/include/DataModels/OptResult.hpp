#pragma once
#ifndef KEVOOPT_DATAMODELS_OPTRESULT_H_INCLUDE
#define KEVOOPT_DATAMODELS_OPTRESULT_H_INCLUDE

#include <vector>

struct OptResult
{
  float fitness;
  std::vector<float> individual;
};

#endif //HIFUEVO_INCLUDE_OPTIMIZER_OPTRESULT_H_
