#include "include/parameters.h"
#include "include/cmaes.h"

#include <cmath>
#include <limits>
#include <ostream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace CMAESpp
{

template<typename T>
Parameters<T>::Parameters()
    : N(-1),
      xstart(0),
      typicalX(0),
      typicalXcase(false),
      rgInitialStds(0),
      rgDiffMinChange(0),
      stopMaxFunEvals(-1),
      facmaxeval(1.0),
      stopMaxIter(-1.0),
      stopTolFun(1e-12),
      stopTolFunHist(1e-13),
      stopTolX(0), // 1e-11*insigma would also be reasonable
      stopTolUpXFactor(1e3),
      lambda(-1),
      mu(-1),
      mucov(-1),
      mueff(-1),
      weights(0),
      damps(-1),
      cs(-1),
      ccumcov(-1),
      ccov(-1),
      facupdateCmode(1),
      weightMode(UNINITIALIZED_WEIGHTS),
      resumefile(""),
      logWarnings(false),
      logStream(std::cerr)
{
  stStopFitness.flg   = false;
  stStopFitness.val   = -std::numeric_limits<T>::max();
  updateCmode.modulo  = -1;
  updateCmode.maxtime = -1;
}

template<typename T>
Parameters<T>::Parameters(const Parameters &parameters) : logStream(std::cerr)
{
  assign(parameters);
}

template<typename T>
Parameters<T>::~Parameters()
{
  if (xstart)
  {
    delete[] xstart;
  }
  if (typicalX)
  {
    delete[] typicalX;
  }
  if (rgInitialStds)
  {
    delete[] rgInitialStds;
  }
  if (rgDiffMinChange)
  {
    delete[] rgDiffMinChange;
  }
  if (weights)
  {
    delete[] weights;
  }
}

template<typename T>
Parameters<T> &Parameters<T>::operator=(const Parameters<T> &parameters)
{
  assign(parameters);
  return *this;
}

/**
 * @param dimension Dimension of the search space \f$N\f$. No default
 *                  available, must be defined here or you have to set the
 *                  member manually.
 * @param inxstart Initial point in search space \f$x_0\f$, default (NULL) is
 *                 \f$(0.5,\ldots,0.5)^T + N(0, initialStdDev^2) \in R^N\f$.
 *                 This must be an array of size \f$N\f$.
 * @param inrgsigma Coordinatewise initial standard deviation of the sample
 *                  distribution (\f$\sigma \cdot \sqrt{C_{ii}} =
 *                  initialStdDev[i]\f$). The expected initial distance
 *                  between initialX and the optimum per coordinate should be
 *                  roughly initialStdDev. The entries should not differ by
 *                  several orders of magnitude. Default (NULL) is
 *                  \f$(0.3,\ldots,0.3)^T \in R^N\f$. This must be an array of
 *                  size \f$N\f$.
 */
template<typename T>
void Parameters<T>::init(int dimension, const T *inxstart, const T *inrgsigma)
{
  if (logWarnings)
  {
    if (!(xstart || inxstart || typicalX))
    {
      logStream << "Warning: initialX undefined. typicalX = 0.5...0.5." << std::endl;
    }
    if (!(rgInitialStds || inrgsigma))
    {
      logStream << "Warning: initialStandardDeviations undefined. 0.3...0.3." << std::endl;
    }
  }

  if (dimension <= 0 && N <= 0)
  {
    throw std::runtime_error("Problem dimension N undefined.");
  }
  else if (dimension > 0)
  {
    N = dimension;
  }

  if (weightMode==UNINITIALIZED_WEIGHTS)
  {
    weightMode = LOG_WEIGHTS;
  }

  diagonalCov = 0; // default is 0, but this might change in future

  if (!xstart)
  {
    xstart = new T[N];
    if (inxstart)
    {
      for (int i = 0; i < N; ++i)
      {
        xstart[i] = inxstart[i];
      }
    }
    else if (typicalX)
    {
      typicalXcase = true;
      for (int i = 0; i < N; ++i)
      {
        xstart[i] = typicalX[i];
      }
    }
    else
    {
      typicalXcase = true;
      for (int i = 0; i < N; i++)
      {
        xstart[i] = 0.5;
      }
    }
  }

  if (!rgInitialStds)
  {
    rgInitialStds = new T[N];
    if (inrgsigma)
    {
      for (int i = 0; i < N; ++i)
      {
        rgInitialStds[i] = inrgsigma[i];
      }
    }
    else
    {
      for (int i = 0; i < N; ++i)
      {
        rgInitialStds[i] = T(0.3);
      }
    }
  }

  supplementDefaults();
}

template<typename T>
void Parameters<T>::assign(const Parameters &p)
{
  N = p.N;
  if (xstart)
  {
    delete[] xstart;
  }
  if (p.xstart)
  {
    xstart = new T[N];
    for (int i = 0; i < N; i++)
    {
      xstart[i] = p.xstart[i];
    }
  }

  if (typicalX)
  {
    delete[] typicalX;
  }
  if (p.typicalX)
  {
    typicalX = new T[N];
    for (int i = 0; i < N; i++)
    {
      typicalX[i] = p.typicalX[i];
    }
  }

  typicalXcase = p.typicalXcase;

  if (rgInitialStds)
  {
    delete[] rgInitialStds;
  }
  if (p.rgInitialStds)
  {
    rgInitialStds = new T[N];
    for (int i = 0; i < N; i++)
    {
      rgInitialStds[i] = p.rgInitialStds[i];
    }
  }

  if (rgDiffMinChange)
  {
    delete[] rgDiffMinChange;
  }
  if (p.rgDiffMinChange)
  {
    rgDiffMinChange = new T[N];
    for (int i = 0; i < N; i++)
    {
      rgDiffMinChange[i] = p.rgDiffMinChange[i];
    }
  }

  stopMaxFunEvals = p.stopMaxFunEvals;
  facmaxeval      = p.facmaxeval;
  stopMaxIter     = p.stopMaxIter;

  stStopFitness.flg = p.stStopFitness.flg;
  stStopFitness.val = p.stStopFitness.val;

  stopTolFun       = p.stopTolFun;
  stopTolFunHist   = p.stopTolFunHist;
  stopTolX         = p.stopTolX;
  stopTolUpXFactor = p.stopTolUpXFactor;

  lambda = p.lambda;
  mu     = p.mu;
  mucov  = p.mucov;
  mueff  = p.mueff;

  if (weights)
  {
    delete[] weights;
  }
  if (p.weights)
  {
    weights = new T[mu];
    for (int i = 0; i < mu; i++)
    {
      weights[i] = p.weights[i];
    }
  }

  damps       = p.damps;
  cs          = p.cs;
  ccumcov     = p.ccumcov;
  ccov        = p.ccov;
  diagonalCov = p.diagonalCov;

  updateCmode.modulo  = p.updateCmode.modulo;
  updateCmode.maxtime = p.updateCmode.maxtime;

  facupdateCmode = p.facupdateCmode;

  weightMode = p.weightMode;

  resumefile = p.resumefile;
}

/**
 * Supplements default parameter values.
 */
template<typename T>
void Parameters<T>::supplementDefaults()
{
  if (lambda < 2)
  {
    lambda = 4 + (int) (3.0*log((double) N));
  }
  if (mu <= 0)
  {
    mu = lambda/2;
  }
  if (!weights)
  {
    setWeights(weightMode);
  }

  if (cs > 0)
  {
    cs *= (mueff + 2.)/(N + mueff + 3.);
  }
  if (cs <= 0 || cs >= 1)
  {
    cs = (mueff + 2.)/(N + mueff + 3.);
  }

  if (ccumcov <= 0 || ccumcov > 1)
  {
    ccumcov = 4./(N + 4);
  }

  if (mucov < 1)
  {
    mucov = mueff;
  }
  T t1 = 2./((N + 1.4142)*(N + 1.4142));
  T t2 = (2.*mueff - 1.)/((N + 2.)*(N + 2.) + mueff);
  t2    = (t2 > 1) ? 1 : t2;
  t2    = (1./mucov)*t1 + (1. - 1./mucov)*t2;
  if (ccov >= 0)
  {
    ccov *= t2;
  }
  if (ccov < 0 || ccov > 1)
  {
    ccov = t2;
  }

  if (diagonalCov < 0)
  {
    diagonalCov = 2 + 100.*N/sqrt((double) lambda);
  }

  if (stopMaxFunEvals <= 0)
  {
    stopMaxFunEvals = facmaxeval*900*(N + 3)*(N + 3);
  }
  else
  {
    stopMaxFunEvals *= facmaxeval;
  }

  if (stopMaxIter <= 0)
  {
    stopMaxIter = ceil((double) (stopMaxFunEvals/lambda));
  }

  if (damps < T(0))
  {
    damps = T(1);
  }
  damps = damps
      *(T(1) + T(2)*std::max(T(0), std::sqrt((mueff - T(1))/(N + T(1))) - T(1)))
      *(T) std::max(T(0.3), T(1) - // modify for short runs
          (T) N/(T(1e-6) + std::min(stopMaxIter, stopMaxFunEvals/lambda)))
      + cs;

  if (updateCmode.modulo < 0)
  {
    updateCmode.modulo = 1./ccov/(double) N/10.;
  }
  updateCmode.modulo *= facupdateCmode;
  if (updateCmode.maxtime < 0)
  {
    updateCmode.maxtime = 0.20;
  } // maximal 20% of CPU-time
}

/**
 * Initializes the offspring weights.
 */
template<typename T>
void Parameters<T>::setWeights(Weights mode)
{
  if (weights)
  {
    delete[] weights;
  }
  weights = new T[mu];
  switch (mode)
  {
    case LINEAR_WEIGHTS:
      for (int i = 0; i < mu; ++i)
      {
        weights[i] = mu - i;
      }
      break;
    case EQUAL_WEIGHTS:
      for (int i = 0; i < mu; ++i)
      {
        weights[i] = 1;
      }
      break;
    case LOG_WEIGHTS:
    default:
      for (int i = 0; i < mu; ++i)
      {
        weights[i] = log(mu + 1.) - log(i + 1.);
      }
      break;
  }

  // normalize weights vector and set mueff
  T        s1 = 0, s2 = 0;
  for (int i  = 0; i < mu; ++i)
  {
    s1 += weights[i];
    s2 += weights[i]*weights[i];
  }
  mueff       = s1*s1/s2;
  for (int i = 0; i < mu; ++i)
  {
    weights[i] /= s1;
  }

  if (mu < 1 || mu > lambda || (mu==lambda && weights[0]==weights[mu - 1]))
  {
    throw std::runtime_error("setWeights(): invalid setting of mu or lambda");
  }
}
}


template class CMAESpp::Parameters<float>  ;
template class CMAESpp::Parameters<double> ;
