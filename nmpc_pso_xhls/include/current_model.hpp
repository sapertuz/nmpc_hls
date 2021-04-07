#ifndef CURRENT_MODEL_HPP
#define CURRENT_MODEL_HPP

#include "system.hpp"

#ifdef INVERTED_PENDULUM_CONFIG
    #include "inverted_pendulum.hpp"
    typedef InvertedPendulum ModelState;
#elif defined SNIFFBOT_CONFIG
    #include "sniffbot.hpp"
    typedef Sniffbot ModelState;
#endif

#endif