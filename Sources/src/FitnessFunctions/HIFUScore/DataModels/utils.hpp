#pragma once
#ifndef _UTILS_H_INCLUDED_
#define _UTILS_H_INCLUDED_

  #ifndef SEP
    #ifdef _WIN32
      #define SEP '\\'
    #else 
      #define SEP '/'
    #endif
  #endif

  #ifndef __INTEL_COMPILER
    #ifdef _WIN32
      #include <malloc.h>
      #define _mm_malloc(a, b) _aligned_malloc((a), (b))
      #define _mm_free(p) _aligned_free(p)
    #else
      #include <stdlib.h>
      #define _mm_malloc(a, b) aligned_alloc((b), (a))  //c++11
      #define _mm_free(p) free(p)
    #endif
  #endif

  #ifndef _ALIGN_
    #define _ALIGN_ 64
  #endif

  #ifdef _WIN32
    #ifndef M_PI
      #define M_PI 3.14159265358979323846
    #endif
  #endif

#endif