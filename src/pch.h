#pragma once

#include "Bc7Mode.h"

#if defined(WIN32)
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdint.h>
#include <string.h>
#include <math.h>

#if defined(OPTION_COUNTERS)
#include <stdio.h>
#include <atomic>
#endif

#if defined(OPTION_AVX512)
#include <immintrin.h> // AVX512
#elif defined(OPTION_AVX2)
#include <immintrin.h> // AVX2
#else
#include <smmintrin.h> // SSE4.1
#endif

#if defined(WIN32)
#include <stdlib.h>
#endif

#if defined(WIN32)
#define INLINED __forceinline
#define NOTINLINED __declspec(noinline)
#else
#define INLINED __attribute__((always_inline))
#define NOTINLINED __attribute__((noinline))
#endif

#define PRINTF(...) printf(__VA_ARGS__); printf("\n");

#if defined(OPTION_AVX512)
constexpr __mmask8 kFullMask8 = 0xFFui8;
constexpr __mmask16 kFullMask16 = 0xFFFFui16;
constexpr __mmask32 kFullMask32 = ~0ui32;
constexpr __mmask64 kFullMask64 = ~0ui64;
#endif
