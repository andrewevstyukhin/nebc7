
//#define OPTION_AVX512
#define OPTION_AVX2
//#define OPTION_FMA
//#define OPTION_PCA
//#define OPTION_COUNTERS
//#define OPTION_SLOWPOKE
#define OPTION_SELFCHECK

//enum { kDenoise = 0, kDenoiseStep = 0 };
//enum { kDenoise = 0, kDenoiseStep = 3 * 3 };
//enum { kDenoise = 1, kDenoiseStep = 0 };
enum { kDenoise = 1, kDenoiseStep = 3 * 3 };

enum { kDenoiseShift = kDenoise ? kDenoise : 1 };

constexpr int kWeightLimit = 1000;
extern int gWeightAlpha, gWeightGreen, gWeightRed, gWeightBlue;
extern int gWeightColor, gWeightColorAlpha;

#if defined(OPTION_AVX512) && (!defined(__AVX512F__) || !defined(__AVX512BW__) || !defined(__AVX512VL__) || defined(OPTION_SLOWPOKE))
#error AVX-512 is required
#endif

#if defined(OPTION_AVX512) && !defined(OPTION_AVX2)
#define OPTION_AVX2
#endif

#if defined(OPTION_AVX2) && (!defined(__AVX2__) || defined(OPTION_SLOWPOKE))
#error AVX2 is required
#endif

#if defined(OPTION_AVX2) && !defined(OPTION_FMA) // Except Via Cores
#define OPTION_FMA
#endif
