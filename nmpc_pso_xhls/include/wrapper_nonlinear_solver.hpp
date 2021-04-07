#ifndef WRAPPER_NSOLVER_HPP
#define WRAPPER_NSOLVER_HPP

#include "nonlinear_solver.hpp"
#include "pso.hpp"
#include "current_model.hpp"

NonlinearSolver * solver = new PSO();

int run_nonlinear_solver(
    _real x_curr[_Nx], 
    _real u_curr[_n_U], 
    int iteration, 
    _real last_best[_Nu*_n_U], 
    _real xref[_Nu*_Nx], 
    _real uref[_n_U], 
    _real xss[_Nx],
    _real uss[_n_U], 
    _real new_best[_Nu*_n_U],
    _real * J
){
    return solver->execute(x_curr, u_curr, iteration, last_best, xref, uref, xss, uss, new_best, J);
}

#endif