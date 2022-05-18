#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#ifdef __SYNTHESIS__
#include "hls_math.h"
typedef half _real_model;
#else
typedef float _real_model;
#endif

#include <time.h>
#include <sys/time.h>

// Definitions
#define N_DOFS 7
#define BILLION 1E9

#if defined(_POSIX_MONOTONIC_CLOCK)
/*  The identifier for the system-wide monotonic clock, which is defined
 *  as a clock whose value cannot be set via clock_settime() and which
 *  cannot have backward clock jumps. */

    #define CLOCK_ID CLOCK_MONOTONIC
#else
    #define CLOCK_ID CLOCK_REALTIME
#endif


#ifdef INVERTED_PENDULUM_CONFIG
#include "hls_inverted_pendulum.hpp"
typedef model_inverted_pendulum<_real_model> model_t;
#elif defined(SNIFFBOT_CONFIG)
#include "hls_sniffbot.hpp"
typedef model_sniffbot<_real_model> model_t;
#endif

#endif
