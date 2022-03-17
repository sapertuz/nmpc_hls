#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#ifdef __SYNTHESIS__
#include "hls_math.h"
typedef half _real_model;
#else
typedef float _real_model;
#endif

#ifdef INVERTED_PENDULUM_CONFIG
#include "hls_inverted_pendulum.hpp"
typedef model_inverted_pendulum<_real_model> model_t;
#elif defined(SNIFFBOT_CONFIG)
#include "hls_sniffbot.hpp"
typedef model_sniffbot<_real_model> model_t;
#endif

#endif
