#ifndef INVERTED_PENDULUM_HPP
#define INVERTED_PENDULUM_HPP

#include "config.hpp"
#include "system.hpp"

class InvertedPendulum : public System {

public:
	InvertedPendulum();
	void model(
        _real state_dot[_Nx],
        _real state[_Nx],
        _real control[_n_U]
    );
    _real nmpc_cost_function(
        _real control_guess[_n_U*_Nu],
        _real xref[_n_U*_Nu], 
        _real uref[_n_U], 
        _real xss[_Nx], 
        _real uss[_n_U]
    );
};

#endif