#ifndef HLS_NONSOLV_HPP
#define HLS_NONSOLV_HPP

/***************************** Include Files *********************************/
#include "config.hpp"
#include "hls_pseudorand.hpp"
#include "hls_system.hpp"
#include "hls_pso.hpp"

/************************** Class Definitions ****************************/

typedef _hw_top_real _pso_real;

int nonlinear_solver_wrapper(
    volatile float *x_curr,//[_Nx], 
    volatile float *u_curr,//[_n_U], 
    int iteration, 
    volatile float *last_best,//[_Nu*_n_U], 
    volatile float *xref,//[_Nu*_Nx], 
    volatile float *uref,//[_n_U], 
    volatile float *xss,//[_Nx],
    volatile float *uss,//[_n_U], 
    
    float *new_best,//[_Nu*_n_U],
    float *J
);

#endif