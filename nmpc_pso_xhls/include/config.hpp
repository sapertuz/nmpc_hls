#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <stdexcept>
    
typedef unsigned short _uchar; //typedef double _real;
typedef float _real; //typedef double _real;
#ifdef __SYNTHESIS__
    #include <hls_math.h>
    typedef half _hw_real;
#else
    typedef float _hw_real;
#endif

// Model configuration variables
//#define INVERTED_PENDULUM_CONFIG
//#define SNIFFBOT_CONFIG

#ifdef INVERTED_PENDULUM_CONFIG
    
    #define _N      20      // Prediction Horizon
    #define _Nu     20      // Control Horizon
    #define _Nx     4       // Number of States
    #define _n_U    1       // Number of Inputs

    const _real  _Ts = 0.1;
    const _real  _u_max[]  = {50.0};
    const _real  _u_min[]  = {-50.0};
    const _real  _du_max[] = {50.0};
    const _real  _uss = 0.0;

    #define _Parametrization 0
    const _real _Lambda = 10;
    const _real _q_param = 10;
    const _real _pmax[] = {1.0};
    const _real _pmin[] = {-1.0};
    
    const _real _state_upper_limits[] = {0.5, 1e10, 1e10, 1e10} ;
    const _real _state_lower_limits[] = {-0.5, -1e10, -1e10, -1e10} ;
    const _real _Q[] = {1e4, 1.0, 1e4, 1.0};
    const _real _Qf[] = {1e3, 1e3, 1e3, 1e3};
    const _real _R[] = {1.0};
    const _real _current_state[] = {0.0, 0.0, 3.1415926536, 0.0};

    #define _Rising_Time 0
    const _real _tr[] =  {0, 0, 0, 0};
#elif defined SNIFFBOT_CONFIG
    #define _N 50
    #define _Nu 50
    #define _Nx 12
    #define _n_U 4
    const _real  _Ts = 0.05;
    const _real  _u_max[] =  {100, 100, 100, 100};
    const _real  _u_min[] =  {-100, -100, -100, -100};
    const _real  _du_max[] =  {20, 20, 20, 20};
    const _real  _uss = 0.0;

    #define  _Parametrization 0
    const _real  _Lambda = 5;
    const _real  _q_param = 2;
    const _real  _pmax[] =  {1, 1, 1, 1};
    const _real  _pmin[] =  {-1, -1, -1, -1};
    
    const _real  _state_upper_limits[] =  {1e10, 1e10, 1e10, 1e10, 1e10, 1e10, 1e10, 1e10, 1e10, 1e10, 1e10, 1e10};
    const _real  _state_lower_limits[] =  {-1e10, -1e10, -1e10, -1e10, -1e10, -1e10, -1e10, -1e10, -1e10, -1e10, -1e10, -1e10};
    const _real  _Q[] =  {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
    const _real  _Qf[] =  {10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0};
    const _real  _R[] =  {0.02, 0.0, 0.0, 0.0};
    const _real  _current_state[] =  {7.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    #define  _Rising_Time 0
    const _real _tr[] =  {10, 10, 10, 0, 0, 8, 0, 0, 0, 0, 0, 0};
#endif

// PSO configuration variables
#define _S 10
#define _KPSO 1
#define _stable_zero 1
#define _maxiter 100
#define _max_v 10
#define _w0 0.9
#define _wf 0.1
#define _c1 2.1
#define _c2 1.0
#define _threshold 1e-2
#define _stop_criteria 0

//#define FAST_SINCOS

#ifdef _WIN32
   //define something for Windows (32-bit and 64-bit, this part is common)
   #ifdef _WIN64
      //define something for Windows (64-bit only)
   #else
      //define something for Windows (32-bit only)
   #endif
#elif __APPLE__

#elif __linux__
    // linux
#elif __unix__ // all unices not caught above
    // Unix
#elif defined(_POSIX_VERSION)
    // POSIX
#else
#   error "Unknown compiler"
#endif

#define BILLION 1E9

#endif
