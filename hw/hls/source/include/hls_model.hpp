#ifndef __MODEL_HPP__
#define __MODEL_HPP__

#ifdef INVERTED_PENDULUM_CONFIG
#include "hls_inverted_pendulum.hpp"
#elif defined(SNIFFBOT_CONFIG)
#include "hls_sniffbot.hpp"
#elif defined(PANDA_CONFIG)
#include "hls_panda.hpp"
#endif

#endif