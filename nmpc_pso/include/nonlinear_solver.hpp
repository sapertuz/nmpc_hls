#ifndef NONLINEAR_SOLVER_HPP
#define NONLINEAR_SOLVER_HPP

#include <iostream>
#include "config.h"
#include "system.hpp"
#include "plot.hpp"

class NonlinearSolver {

protected: 
	System * controled_system;

public:
	virtual ~NonlinearSolver(){};
	virtual int execute(_real * u_curr, int iteration, _real * last_best, _real * xref, _real * uref, _real * xss, _real * uss, _real * J, Plot * graph) = 0;
};

#endif
