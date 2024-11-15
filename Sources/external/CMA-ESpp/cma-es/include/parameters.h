#pragma once

#include <iostream>
#include <string>

namespace CMAESpp
{
template<typename T>
class CMAES;

/**
 * @class Parameters
 * Holds all parameters that can be adjusted by the user.
 */
template<typename T>
class Parameters
{
  friend class CMAES<T>;

 public:

  /* Input parameters. */
  //! Problem dimension, must stay constant.
  int N;
  //! Initial search space vector.
  T *xstart;
  //! A typical value for a search space vector.
  T *typicalX;
  //! Indicates that the typical x is the initial point.
  bool typicalXcase;
  //! Initial standard deviations.
  T *rgInitialStds;
  T *rgDiffMinChange;

  /* Termination parameters. */
  //! Maximal number of objective function evaluations.
  T stopMaxFunEvals;
  T facmaxeval;
  //! Maximal number of iterations.
  T stopMaxIter;
  //! Minimal fitness value. Only activated if flg is true.
  struct
  {
    bool flg;
    T    val;
  } stStopFitness;
  //! Minimal value difference.
  T stopTolFun;
  //! Minimal history value difference.
  T stopTolFunHist;
  //! Minimal search space step size.
  T stopTolX;
  //! Defines the maximal condition number.
  T stopTolUpXFactor;

  /* internal evolution strategy parameters */
  /**
   * Population size. Number of samples per iteration, at least two,
   * generally > 4.
   */
  int lambda;
  /**
   * Number of individuals used to recompute the mean.
   */
  int mu;
  T   mucov;
  /**
   * Variance effective selection mass, should be lambda/4.
   */
  T   mueff;
  /**
   * Weights used to recombinate the mean sum up to one.
   */
  T *weights;
  /**
   * Damping parameter for step-size adaption, d = inifinity or 0 means adaption
   * is turned off, usually close to one.
   */
  T damps;
  /**
   * cs^-1 (approx. n/3) is the backward time horizon for the evolution path
   * ps and larger than one.
   */
  T cs;
  T ccumcov;
  /**
   * ccov^-1 (approx. n/4) is the backward time horizon for the evolution path
   * pc and larger than one.
   */
  T ccov;
  T diagonalCov;
  struct
  {
    T modulo;
    T maxtime;
  } updateCmode;
  T facupdateCmode;

  /**
   * Determines the method used to initialize the weights.
   */
  enum Weights
  {
    UNINITIALIZED_WEIGHTS, LINEAR_WEIGHTS, EQUAL_WEIGHTS, LOG_WEIGHTS
  } weightMode;

  //! File that contains an optimization state that should be resumed.
  std::string resumefile;

  //! Set to true to activate logging warnings.
  bool logWarnings;
  //! Output stream that is used to log warnings, usually std::cerr.
  std::ostream &logStream;

  Parameters();

  Parameters(const Parameters<T> &parameters);

  ~Parameters();

  Parameters<T> &operator=(const Parameters<T> &parameters);

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
  void init(int dimension = 0, const T *inxstart = 0, const T *inrgsigma = 0);

 private:
  void assign(const Parameters<T> &p);

  /**
   * Supplements default parameter values.
   */
  void supplementDefaults();

  /**
   * Initializes the offspring weights.
   */
  void setWeights(Weights mode);
};
}