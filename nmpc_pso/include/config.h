#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <stdexcept>

//typedef double _real;
typedef float _real;

//#define DEBUG
//#define PRINT_TO_TERMINAL
//#define PLOT
//#define SAVE_PLOT

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define _N  20
#define _Nu 20
#define _Nx 4

//#define GNUPLOT

//#define FAST_SINCOS

#ifdef _WIN32
   //define something for Windows (32-bit and 64-bit, this part is common)
	#define GNUPLOT_TERM "qt"
   #ifdef _WIN64
      //define something for Windows (64-bit only)
   #else
      //define something for Windows (32-bit only)
   #endif
#elif __APPLE__
	#define GNUPLOT_TERM "aqua"
    // #include "TargetConditionals.h"
    // #if TARGET_IPHONE_SIMULATOR
    //      // iOS Simulator
    // #elif TARGET_OS_IPHONE
    //     // iOS device
    // #elif TARGET_OS_MAC
    //     // Other kinds of Mac OS
    // #else
    // #   error "Unknown Apple platform"
    // #endif
#elif __linux__
	#define GNUPLOT_TERM "qt"
    // linux
#elif __unix__ // all unices not caught above
    // Unix
	#define GNUPLOT_TERM "qt"
#elif defined(_POSIX_VERSION)
    // POSIX
	#define GNUPLOT_TERM "qt"
#else
#   error "Unknown compiler"
#endif

#define BILLION 1E9

#endif
