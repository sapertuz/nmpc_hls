#ifndef SNIFFBOT_HPP
#define SNIFFBOT_HPP

#include "config.hpp"
#include "system.hpp"

class Sniffbot : public System {

public:
	Sniffbot();
	void model(
        _real state_dot[_Nx],
        _real state[_Nx],
        _real control[_n_U]
    );
    _real nmpc_cost_function(
        _real control_guess[_n_U*_Nu],
        _real xref[_Nx*_Nu],
        _real uref[_n_U], 
        _real xss[_Nx], 
        _real uss[_n_U]
    );
private:
	_real value_map(
		_real input, 
		_real in_min, 
		_real in_max, 
		_real out_min, 
		_real out_max
	);
};

#endif