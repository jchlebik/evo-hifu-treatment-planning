/**
 * @file      kWaveDiffusionSolver2D.cpp
 *
 * @author    Filip Kuklis \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            ikuklis@fit.vutbr.cz
 *
 * @brief     The implementation file containing the main class of the project responsible for the simulation.
 *
 * @version   kWaveDiffusion
 *
 * @date      05 November      2018, 09:45 (created) \n
 *            12 November      2018, 15:00 (revised)
 *
 * @copyright Copyright (C) 2018 Filip Kuklis and Bradley Treeby.
 *
 * This file is part of the C++ extension of the [k-Wave Toolbox](http://www.k-wave.org).
 *
 * This file is part of the k-Wave. k-Wave is free software: you can redistribute it and/or modify it under the terms
 * of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with k-Wave.
 * If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
 */



// Linux build
#ifdef __linux__
#include <sys/resource.h>
#endif

// Windows build
#ifdef _WIN64
#define _USE_MATH_DEFINES
  #include <Windows.h>
  #include <Psapi.h>
  #pragma comment(lib, "Psapi.lib")
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include "kWaveDiffusionSolver.h"
//#include "../Logger/Logger.h"
//#include <immintrin.h>
#include <cmath>
#include <ctime>
#include <limits>
#include <fftw3.h>

#include <stdlib.h>

#ifndef __INTEL_COMPILER
  #define _mm_malloc(a, b) aligned_alloc((b), (a))
  #define _mm_free(p) free(p)
  #define __assume_aligned(p, n); 
#endif

#define _ALIGN_ 64

using std::ios;

//--------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------- Constants -----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor of the class.
 */
kWaveDiffusionSolver2D::kWaveDiffusionSolver2D(int NX, int NY, const float* Temperature, const float* mediumDensity, const float* mediumSpecificHeat, float* objQ, const float* objCem43, const float* mediumThermalConductivity, float* output1, float* output2)  :
            mMatlabObjTemperature(Temperature), mMatlabObjQ(objQ), mMatlabObjCem43(objCem43), mMatlabMediumDensity(mediumDensity), mMatlabMediumSpecificHeat(mediumSpecificHeat), mMatlabMediumThermalConductivity(mediumThermalConductivity), mMatlabOutputPtr1(output1), mMatlabOutputPtr2(output2),
            mQScaleFac(0),  /* TODO mQScaleFac is not used now*/ mDdiffusionCoeffRef(0), mNt(0), mUseQterm(0), mDt(0)
           //mParameters(CommandLineParameters::getInstance())
{
  mNx = NX;
  mNxr = (NX + 1) / 2;
  mNy = NY;
  //auto a = omp_get_max_threads();
  omp_set_num_threads(omp_get_max_threads());

  try
  {
      allocateMemory();
  }
  catch (std::bad_alloc &ba)
  {
      // Logger::log(Logger::LogLevel::kBasic, kOutFmtFailed);
      // Logger::log(Logger::LogLevel::kBasic, kOutFmtLastSeparator);
      // Logger::errorAndTerminate(Logger::wordWrapString(ba.what(),kErrFmtPathDelimiters, 9));
  }

  preProcessing();

  InitializeFftwPlans();

  // does not load data from disk, only copy Matlab created matrices(those which need to be changed) to allocated matrices
  // not necessery if you are responsible for matrices and do not need to copy them (you can change them)
  loadInputData();

  createDiffusionTerms();
}// end of kWaveDiffusionSolver2D
//----------------------------------------------------------------------------------------------------------------------


/**
 * Destructor of the class.
 */
kWaveDiffusionSolver2D::~kWaveDiffusionSolver2D()
{
    freeMemory();
}// end of kWaveDiffusionSolver2D
//----------------------------------------------------------------------------------------------------------------------

/**
 * The method allocates the matrix container, creates all matrices and creates all output streams
 * (however not allocating memory).
 */
void kWaveDiffusionSolver2D::allocateMemory()
{
  mTempRealArray1 =               (float *) _mm_malloc(mNx  * mNy * sizeof(float),               _ALIGN_);
  mTempRealArray2 =               (float *) _mm_malloc(mNx  * mNy * sizeof(float),               _ALIGN_);
  mQTerm          =               (float *) _mm_malloc(mNx  * mNy * sizeof(float),               _ALIGN_);
  mPTerm          =               (float *) _mm_malloc(mNx  * mNy * sizeof(float),               _ALIGN_);
  mDiffusionP1    =               (float *) _mm_malloc(mNx  * mNy * sizeof(float),               _ALIGN_);
  mDiffusionP2    =               (float *) _mm_malloc(mNx  * mNy * sizeof(float),               _ALIGN_);
  mTemperature    =               (float *) _mm_malloc(mNx  * mNy * sizeof(float),               _ALIGN_);
  mCem43          =               (float *) _mm_malloc(mNx  * mNy * sizeof(float),               _ALIGN_);
  mKx_vec         =               (float *) _mm_malloc(mNxr       * sizeof(float),               _ALIGN_);
  mKy_vec         =               (float *) _mm_malloc(       mNy * sizeof(float),               _ALIGN_);
  mKappa          =               (float *) _mm_malloc(mNxr * mNy * sizeof(float),               _ALIGN_);
  mK2             =               (float *) _mm_malloc(mNxr * mNy * sizeof(float),               _ALIGN_);
  mTempFftwArray1 = (std::complex<float> *) _mm_malloc(mNxr * mNy * sizeof(std::complex<float>), _ALIGN_);
  mTempFftwArray2 = (std::complex<float> *) _mm_malloc(mNxr * mNy * sizeof(std::complex<float>), _ALIGN_);
  mDerivX         = (std::complex<float> *) _mm_malloc(mNxr * mNy * sizeof(std::complex<float>), _ALIGN_);
  mDerivY         = (std::complex<float> *) _mm_malloc(mNxr * mNy * sizeof(std::complex<float>), _ALIGN_);
}// end of allocateMemory
//----------------------------------------------------------------------------------------------------------------------


/**
 * The method frees all memory allocated by the class.
 */
void kWaveDiffusionSolver2D::freeMemory()
{
  _mm_free(mTempRealArray1);
  _mm_free(mTempRealArray2);
  _mm_free(mQTerm);
  _mm_free(mPTerm);
  _mm_free(mDiffusionP1);
  _mm_free(mDiffusionP2);
  _mm_free(mTemperature);
  _mm_free(mCem43);
  _mm_free(mKappa);
  _mm_free(mK2); //-
  _mm_free(mKx_vec);
  _mm_free(mKy_vec);
  _mm_free(mTempFftwArray1); //-
  _mm_free(mTempFftwArray2); //-
  _mm_free(mDerivX); //-
  _mm_free(mDerivY); //-
}// end of freeMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * Get peak memory usage.
 */
size_t kWaveDiffusionSolver2D::getMemoryUsage() const
{
    // Linux build
#ifdef __linux__
    struct rusage mem_usage;
    getrusage(RUSAGE_SELF, &mem_usage);

    return mem_usage.ru_maxrss >> 10;
#endif

    // Windows build
#ifdef _WIN64
    HANDLE hProcess;
    PROCESS_MEMORY_COUNTERS pmc;

    hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ,
                           FALSE,
                           GetCurrentProcessId());

    GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc));
    CloseHandle(hProcess);

    return pmc.PeakWorkingSetSize >> 20;
#endif
}// end of getMemoryUsage
//----------------------------------------------------------------------------------------------------------------------

/**
 * Get release code version.
 */
std::string kWaveDiffusionSolver2D::getCodeName() const
{
  return "NYI";
}// end of getCodeName
//----------------------------------------------------------------------------------------------------------------------


/**
 * Print full code name and the license.
 */
void kWaveDiffusionSolver2D::printFullCodeNameAndLicense() const
{

}// end of printFullCodeNameAndLicense
//----------------------------------------------------------------------------------------------------------------------

/**
 * Copy data from matlab to local allocated memory (mex input array should not be changed)
 */
void kWaveDiffusionSolver2D::loadInputData()
{
  const float* pMatlabObjTemp  = mMatlabObjTemperature;
  const float* pMatlabObjCem43 = mMatlabObjCem43;
  float*       pTemperature    = mTemperature;
  float*       pCem43          = mCem43;

  #pragma omp parallel for
  for(size_t i = 0; i < mNx*mNy; i++)
  {
      pTemperature[i] = pMatlabObjTemp[i];
      pCem43[i]       = pMatlabObjCem43[i];
  }
}// end of loadInputData
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute diffusion terms P1 and P2
 */
void kWaveDiffusionSolver2D::createDiffusionTerms()
{
    const float* pMatlabMediumDensity             = mMatlabMediumDensity;
    const float* pMatlabMediumSpecificHeat        = mMatlabMediumSpecificHeat;
    const float* pMatlabMediumThermalConductivity = mMatlabMediumThermalConductivity;

    float*       pDiffusionP1    = mDiffusionP1;
    float*       pDiffusionP2    = mDiffusionP2;

    #pragma omp parallel for simd aligned(pDiffusionP1, pMatlabMediumDensity, pMatlabMediumSpecificHeat:_ALIGN_)
    for(size_t i = 0; i < mNx*mNy; i++)
    {
        pDiffusionP1[i] = 1.0f / (pMatlabMediumDensity[i] * pMatlabMediumSpecificHeat[i]);
    }

    #pragma omp parallel for simd aligned(pDiffusionP2, pMatlabMediumThermalConductivity:_ALIGN_)
    for(size_t i = 0; i < mNx*mNy; i++)
    {
        pDiffusionP2[i] = pMatlabMediumThermalConductivity[i];
    }

    #pragma omp parallel for reduction(max:mDdiffusionCoeffRef)
    for(size_t i = 0; i < mNx*mNy; i++)
    {
        mDdiffusionCoeffRef = std::max(mDdiffusionCoeffRef, (pDiffusionP1[i] * pDiffusionP2[i]));
    }
}// end of createDiffusionTerms
//----------------------------------------------------------------------------------------------------------------------

/**
* This method computes k-wave diffusion 2D simulation.
 */
void kWaveDiffusionSolver2D::takeTimeStep(int nt, bool useQTerm, float dt)
{
  mNt           = nt;
  mUseQterm     = useQTerm;
  mDt           = dt;

  makeKappa();

  makeDerivMatrix();

  //TODO which qterm
  if(mUseQterm)
  {
    createQTermDiff();
  }

  computeMainLoop();

}// end of compute()
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ Protected methods -------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Initialize FFTW plans.
 */
void kWaveDiffusionSolver2D::InitializeFftwPlans()
{
  fftwf_init_threads();

  fftwf_plan_with_nthreads(omp_get_max_threads());
  mFftwPlanForward = fftwf_plan_dft_r2c_2d(mNx,
                                           mNy,
                                           mTempRealArray1,
                                           reinterpret_cast<fftwf_complex*>(mTempFftwArray1),
                                          FFTW_MEASURE);


  mFftwPlanBackward = fftwf_plan_dft_c2r_2d(mNx,
                                            mNy,
                                            reinterpret_cast<fftwf_complex*>(mTempFftwArray1),
                                            mTempRealArray1,
                                            FFTW_MEASURE);
}// end of InitializeFftwPlans
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute pre-processing phase.
 */
void kWaveDiffusionSolver2D::preProcessing()
{
  float dx = 0.0002; //TODO
  float dy = 0.0002; //TODO

  makeKVec<true>(mKx_vec, dx, mNx);
  makeKVec<false>(mKy_vec, dy, mNy);

  makeK2();
}// end of preProcessing
//----------------------------------------------------------------------------------------------------------------------

template<bool reduced>
void kWaveDiffusionSolver2D::makeKVec(float* k_vec, float d, size_t N)
{
  const float NF = (float) N;

  const float nx = -(((NF - (N % 2 == 0 ? 0.0f : 1.0f)) / 2.0f) / NF);
  const float step = ((-(nx * 2.0f)) / (NF - (N % 2 == 0 ? 0.0f : 1.0f)));

  const float mPi = ((M_PI * 2.0f) / d);

  const size_t half = (N+1) / 2;

  const float mPiNx   = mPi*nx;
  const float mPiStep = mPi*step;

  #pragma omp parallel for simd aligned(k_vec:_ALIGN_)
  for (size_t i = (half - 1); i < N; i++)
  {
    //k_vec[i - (half - 1)] = mPi * (nx + (i * step));
    k_vec[i - (half - 1)] = mPiNx + mPiStep * i;
  }

  if(!reduced)
  {
    #pragma omp parallel for simd aligned(k_vec:_ALIGN_)
    for (size_t i = 0; i < (half-1); i++)
    {
      //k_vec[i + (half)] = mPi * (nx + (i * step));
      k_vec[i + (half)] = mPiNx + mPiStep * i;
    }
  }

  // matlab comment
  //% force middle value to be zero in case 1/Nx is a recurring
  //            % number and the series doesn't give exactly zero
  //            nx(floor(Nx/2) + 1) = 0;
  // matlab comment
  // THIS VECTOR IS SHIFTED THUS MIDDLE VALUE IS AT 0 index
  k_vec[0] = 0.0f;
}

void kWaveDiffusionSolver2D::makeDerivMatrix()
{

  const float*    pKappa    = mKappa;
  const float*    pKx_vec   = mKx_vec;

  #pragma omp parallel for
  for(size_t y = 0; y < mNy; y++)
  {
    const float kyVec = mKy_vec[y];
    const int   indexBase = y*mNxr;

    #pragma omp simd aligned(pKappa, pKx_vec : _ALIGN_)
    for(size_t x = 0; x < mNxr; x++)
    {
      const float sqrtKappa = sqrtf(pKappa[indexBase + x]);
      mDerivX[indexBase + x] = {0.0f, sqrtKappa * pKx_vec[x]};
      mDerivY[indexBase + x] = {0.0f, sqrtKappa * kyVec};
    }
  }
}

void kWaveDiffusionSolver2D::makeK2()
{
  const float*    pKx_vec   = mKx_vec;
  const float*    pKy_vec   = mKy_vec;

  #pragma omp parallel for
  for(size_t y = 0; y < mNy; y++)
  {
    const int indexBase = y*mNxr;
    #pragma omp simd aligned(pKx_vec, pKy_vec : _ALIGN_)
    for(size_t x = 0; x < mNxr; x++)
    {
      mK2[indexBase + x] = pKx_vec[x] * pKx_vec[x] + pKy_vec[y] * pKy_vec[y];
    }
  }
}

void kWaveDiffusionSolver2D::makeKappa()
{
  float perfusionCoeffRef = 0.0f;                  //TODO - perfusion is not used now
  const float diffusionCoeffRef = mDdiffusionCoeffRef;

  float dTdiffusionCoeffRef = mDt * diffusionCoeffRef;

  const float* pK2 = mK2;

  if(perfusionCoeffRef == 0.0f)
  {
    #pragma omp parallel for simd aligned(pK2 :_ALIGN_)
    for (size_t i = 0; i < mNxr*mNy; i++)
    {
      mKappa[i] = dTdiffusionCoeffRef * pK2[i];
      mKappa[i] = (pK2[i] == 0.0f) ? 1.0f : ((1.0f - expf(-mKappa[i]) ) / mKappa[i]);
    }
  }
  else{} //TODO
}




/*  // SHOULD work, but it is NOT TESTED
void kWaveDiffusionSolver2D::createPTerm()
{
  //TODO is abb temp scalar?
  const float ambTemp  = mMatlabObjBlAmbTemp[0];
  const float perfCoef = mMatlabObjPerfcoef[0];
  //TODO

  const float          numOfElem = mNx*mNy;
  std::complex<float>* pTempFftwArray = mTempFftwArray1;
  const float*         pMatlabKappa = mMatlabKappa;
  const float*         pTemperature = mTemperature;
  float*               pPTerm = mPTerm;

#pragma omp parallel for simd
  for(size_t i = 0; i < mNy*mNx; i++)
  {
    pPTerm[i] = pTemperature[i] - ambTemp;
  }


  fftw_execute_dft_r2c(mFftwPlanForward,
                       mPTerm,
                       reinterpret_cast<fftw_complex*>(mTempFftwArray1));

  #pragma omp parallel for simd
  for(size_t i = 0; i < mNy*mNxr; i++)
  {
    // TODO -1.0 * perfCoef ?? look at matlab code
    pTempFftwArray[i] *= -1.0 * perfCoef * pMatlabKappa[i];
    pTempFftwArray[i] /= numOfElem;
  }


  fftw_execute_dft_c2r(mFftwPlanBackward,
                       reinterpret_cast<fftw_complex*>(mTempFftwArray1),
                       mPTerm);
}
*/

/**
 * Compute QTerm using scaling
 */
void kWaveDiffusionSolver2D::createQTermScale()
{
  fftwf_execute_dft_r2c(mFftwPlanForward,
                       mMatlabObjQ,
                       reinterpret_cast<fftwf_complex*>(mTempFftwArray1));

  const float          qScaleFac     = mQScaleFac;
  const float          numOfElemInv  = 1.0f/(mNx*mNy);
  std::complex<float>* tempFftwArray = mTempFftwArray1;
  const float*         pKappa        = mKappa;

  const float          weightFactor  = qScaleFac * numOfElemInv;

  #pragma omp parallel for simd aligned(tempFftwArray, pKappa:_ALIGN_)
  for(size_t i = 0; i < mNy*mNxr; i++)
  {
      tempFftwArray[i] = tempFftwArray[i] * pKappa[i] * weightFactor;
  }

  fftwf_execute_dft_c2r(mFftwPlanBackward,
                       reinterpret_cast<fftwf_complex*>(mTempFftwArray1),
                       mQTerm);
}// end of createQTermScale
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute QTerm using diffusion
 */
void kWaveDiffusionSolver2D::createQTermDiff()
{
  fftwf_execute_dft_r2c(mFftwPlanForward,
                        mMatlabObjQ,
                        reinterpret_cast<fftwf_complex*>(mTempFftwArray1));

  const float          qScaleFac     = mQScaleFac;
  const float          numOfElemInv  = 1.0f/(mNx*mNy);
  std::complex<float>* tempFftwArray = mTempFftwArray1;
  const float*         pKappa        = mKappa;
  float*               pQTerm        = mQTerm;
  const float*         pDiffP1       = mDiffusionP1;


  #pragma omp parallel for simd aligned(tempFftwArray, pKappa:_ALIGN_)
  for(size_t i = 0; i < mNy*mNxr; i++)
  {
    tempFftwArray[i] = tempFftwArray[i] * pKappa[i];
  }

  fftwf_execute_dft_c2r(mFftwPlanBackward,
                        reinterpret_cast<fftwf_complex*>(mTempFftwArray1),
                        mQTerm);

  #pragma omp parallel for simd aligned(pQTerm, pDiffP1:_ALIGN_)
  for(size_t i = 0; i < mNy*mNx; i++)
  {
    pQTerm[i] = pQTerm[i] * pDiffP1[i] * numOfElemInv;
  }
}// end of createQTermDiff
//----------------------------------------------------------------------------------------------------------------------


/**
 * Compute the main time loop of kWaveDiffusionSolver2D.
 */
void kWaveDiffusionSolver2D::computeMainLoop()
{
  const float dt60      = mDt / 60.0f;
  const float numOfElem = (float)mNx*mNy;

  const size_t    iNumOfElem  = mNx*mNy;
  const size_t    iRNumOfElem = mNxr*mNy;

  std::complex<float>* tempFftwArray  = mTempFftwArray1;
  const float*         pQTerm         = mQTerm;
  float*               pTemperature   = mTemperature;
  float*               pTempRealArray = mTempRealArray1;
  float*               pCem43         = mCem43;

  const int            useQT          = mUseQterm;


  //double t1 = omp_get_wtime();

  for(size_t n = 0; n < mNt; n++)
  {
    computeDTermHeterogenous();

    if(useQT)
    {
      #pragma omp parallel for simd aligned(pTemperature, pTempRealArray, pQTerm:_ALIGN_)
      for (size_t i = 0; i < iNumOfElem; i++)
      {
        //float qTerm = (useQT == 1 ? pQTerm[i] : 0.0f);
        pTemperature[i] += mDt * ((pTempRealArray[i]) + /*mPTerm[i] */  pQTerm[i]);
      }
    }
    else
    {
      #pragma omp parallel for simd aligned(pTemperature, pTempRealArray:_ALIGN_)
      for (size_t i = 0; i < iNumOfElem; i++)
      {
        //float qTerm = (useQT == 1 ? pQTerm[i] : 0.0f);
        pTemperature[i] += mDt * ((pTempRealArray[i]));
      }
    }
    #pragma omp parallel for simd aligned(pTemperature, pCem43:_ALIGN_)
    for (size_t i = 0; i < iNumOfElem; i++)
    {
        float half = (pTemperature[i] >= 43.0f ? 0.5f : 0.0f);
        float quart = ((pTemperature[i] >= 37.0f and pTemperature[i] < 43.0f) ? 0.25f : 0.0f);
        float part = std::pow((half + quart),  (43.0f - pTemperature[i]));
        pCem43[i] = pCem43[i] + dt60 * part;
    }
  }

  float* pOutPtr1 = mMatlabOutputPtr1;
  float* pOutPtr2 = mMatlabOutputPtr2;

  #pragma omp parallel for// outPtrs come from outside and are not alligned
  for(size_t i = 0; i < mNx*mNy; i++)
  {
    pOutPtr1[i] = pCem43[i];
    pOutPtr2[i] = pTemperature[i];
  }

}// end of computeMainLoop
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute DTerm for homogenous simulation
 */
void kWaveDiffusionSolver2D::computeDTermHomogenous()
{
  const size_t    iNumOfElem   = mNx*mNy;
  const size_t    iRNumOfElem  = mNxr*mNy;
  const float     numOfElemInv = 1.0f/(mNx*mNy);

  const float          diffP = mDiffusionP1[0] * mDiffusionP1[0];

  const float*         pKappa       = mKappa;
  const float*         pObjK2       = mK2;

  std::complex<float>* pTempFftwArray = mTempFftwArray1;

  fftwf_execute_dft_r2c(mFftwPlanForward,
                       mTemperature,
                       reinterpret_cast<fftwf_complex *>(mTempFftwArray1));

  #pragma omp parallel for simd aligned(pTempFftwArray, pKappa, pObjK2:_ALIGN_)
  for (size_t i = 0; i < iRNumOfElem; i++)
  {
    pTempFftwArray[i] *= pKappa[i] * -pObjK2[i] * diffP * numOfElemInv;
  }

  /// now mTempRealArray1 is dTerm array
  fftwf_execute_dft_c2r(mFftwPlanBackward,
                       reinterpret_cast<fftwf_complex *>(mTempFftwArray1),
                       mTempRealArray1);
}// end of computeDTermHomogenous
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute DTerm for heterogenous simulation
 */
void kWaveDiffusionSolver2D::computeDTermHeterogenous()
{
  const float numOfElemInv = 1.0f/(mNx*mNy);

  std::complex<float>* pTempFftwArray1 = mTempFftwArray1;
  std::complex<float>* pTempFftwArray2 = mTempFftwArray2;

  float*               pTempRealArray1 = mTempRealArray1;
  float*               pTempRealArray2 = mTempRealArray2;


  const std::complex<float>* pDerivX = mDerivX;
  const std::complex<float>* pDerivY = mDerivY;

  // check if diffusion is matrix
  const float*         diffP1 = mDiffusionP1;
  const float*         diffP2 = mDiffusionP2;



  fftwf_execute_dft_r2c(mFftwPlanForward,
                       mTemperature,
                       reinterpret_cast<fftwf_complex*>(mTempFftwArray1));

  #pragma omp parallel for simd aligned(pTempFftwArray2, pTempFftwArray1, pDerivY, pDerivX:_ALIGN_)
  for(size_t i = 0; i < mNy*mNxr; i++)
  {
    /// array1 cant be first
    pTempFftwArray2[i]  = pTempFftwArray1[i] * pDerivY[i] * numOfElemInv;
    pTempFftwArray1[i] *= pDerivX[i] * numOfElemInv;
  }


  fftwf_execute_dft_c2r(mFftwPlanBackward,
                       reinterpret_cast<fftwf_complex*>(mTempFftwArray1),
                       pTempRealArray1);

  fftwf_execute_dft_c2r(mFftwPlanBackward,
                       reinterpret_cast<fftwf_complex*>(mTempFftwArray2),
                       pTempRealArray2);

  #pragma omp parallel for simd aligned(pTempRealArray1, pTempRealArray2, diffP2:_ALIGN_)
  for(size_t i = 0; i < mNy*mNx; i++)
  {
    pTempRealArray1[i] *= diffP2[i];
    pTempRealArray2[i] *= diffP2[i];
  }

  fftwf_execute_dft_r2c(mFftwPlanForward,
                       pTempRealArray1,
                       reinterpret_cast<fftwf_complex*>(mTempFftwArray1));

  fftwf_execute_dft_r2c(mFftwPlanForward,
                       pTempRealArray2,
                       reinterpret_cast<fftwf_complex*>(mTempFftwArray2));

  #pragma omp parallel for simd aligned(pTempFftwArray1, pTempFftwArray2, pDerivX, pDerivY:_ALIGN_)
  for(size_t i = 0; i < mNy*mNxr; i++)
  {
    pTempFftwArray1[i] *= pDerivX[i] * numOfElemInv;
    pTempFftwArray2[i] *= pDerivY[i] * numOfElemInv;
  }

  fftwf_execute_dft_c2r(mFftwPlanBackward,
                       reinterpret_cast<fftwf_complex*>(mTempFftwArray1),
                       pTempRealArray1);

  fftwf_execute_dft_c2r(mFftwPlanBackward,
                       reinterpret_cast<fftwf_complex*>(mTempFftwArray2),
                       pTempRealArray2);

  #pragma omp parallel for simd aligned(pTempRealArray1, pTempRealArray2, diffP1:_ALIGN_)
  for(size_t i = 0; i < mNy*mNx; i++)
  {
    pTempRealArray1[i] = (pTempRealArray1[i] + pTempRealArray2[i]) * diffP1[i];
  }
}// end of computeDTermHeterogenouss
//----------------------------------------------------------------------------------------------------------------------


/**
 * Post processing the quantities, closing the output streams and storing the sensor mask.
 */
void kWaveDiffusionSolver2D::postProcessing()
{

}// end of postProcessing
//----------------------------------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Private methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

