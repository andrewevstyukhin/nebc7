
//#define OPTION_AVX512
#define OPTION_AVX2
//#define OPTION_FMA
//#define OPTION_PCA
//#define OPTION_COUNTERS
//#define OPTION_LINEAR
#define OPTION_SELFCHECK

#if defined(OPTION_LINEAR)

// Linear RGB
enum { kAlpha = 3, kGreen = 1, kRed = 1, kBlue = 1 };

#else

// http://www.brucelindbloom.com/index.html?WorkingSpaceInfo.html sRGB
enum { kAlpha = 1000, kGreen = 715, kRed = 213, kBlue = 72 };

#endif

enum { kColor = kGreen + kRed + kBlue };

#if defined(OPTION_AVX512) && !defined(OPTION_AVX2)
#define OPTION_AVX2
#endif

#if defined(OPTION_AVX2) && !defined(OPTION_FMA) // Except Via Cores
#define OPTION_FMA
#endif

#if defined(OPTION_AVX2) && !defined(__AVX2__)
#error AVX2 is required
#endif

#if !defined(__AVX__)
#error AVX is required
#endif