#ifndef NONLINEAR_SOLVER_HPP
#define NONLINEAR_SOLVER_HPP

#include <iostream>
#include "config.hpp"
#include "system.hpp"

#include "current_model.hpp"

class NonlinearSolver {

protected: 
	System * controled_system = new ModelState();

public:
	virtual ~NonlinearSolver(){};
	virtual int execute(
		_real * x_curr, 
		_real * u_curr, 
		int iteration, 
		_real * last_best,
		_real * xref,
		_real * uref,
		_real * xss,
		_real * uss,
		_real * new_best,
		_real * J
	) = 0;

};

#endif
