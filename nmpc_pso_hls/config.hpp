#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#ifdef INVERTED_PENDULUM_CONFIG
    typedef model_inverted_pendulum<_real> model_t;
#elif defined(SNIFFBOT_CONFIG)
    typedef model_sniffbot<_real> model_t;
#endif

#endif