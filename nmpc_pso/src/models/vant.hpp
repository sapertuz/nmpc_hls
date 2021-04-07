#ifndef VANT_HPP
#define VANT_HPP

#include <string.h>
#include "config.h"
#include "system.hpp"
#include <fstream>

class Vant : public System {

public:
	Vant(std::string config_file);
	void model(_real * state_dot, _real * state, _real * control);
	_real nmpc_cost_function(_real * control_guess, _real * xref, _real * uref, _real * xss, _real * uss);

    void control_from_parameters(_real * parameters, _real * previous_control, _real * control, int horizon);
};

#endif