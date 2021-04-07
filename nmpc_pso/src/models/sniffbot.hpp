#ifndef SNIFFBOT_HPP
#define SNIFFBOT_HPP

#include <string.h>
#include "config.h"
#include "system.hpp"
#include <fstream>

class Sniffbot : public System {

public:
	Sniffbot(std::string config_file);
	void model(_real * state_dot, _real * state, _real * control);
	void control_from_parameters(_real * parameters, _real * previous_control, _real * control, int horizon);
	_real nmpc_cost_function(_real * control_guess, _real * xref, _real * uref, _real * xss, _real * uss);

private:
	_real value_map(_real input, _real in_min, _real in_max, _real out_min, _real out_max);
};

#endif