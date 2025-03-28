#include "hls_nonlinear_solver.hpp"
// #include "hls_pseudorand.hpp"
// #include "hls_system.hpp"
#include "hls_pso.hpp"

/************************** Variable Definitions *****************************/

// Model Core Generator
// top_model_t my_model;
// top_model_t *my_model_ptr = &my_model;
// System Core Generator
// typedef System<_hw_top_real,_Nh, _Nx, _n_U, _Nu> _system_t; 
// _system_t _hw_system(
//     _u_max, 
//     _u_min, 
//     _du_max,
//     _controlled_state,
//     _state_upper_limits, 
//     _state_lower_limits, 
//     _Q, 
//     _Qf, 
//     _R, 
//     _uss, 
//     _Ts
//     // ,
//     // my_model
// );
// Pseudo Random Core Generator
// typedef pseudoRand_gen<_hw_top_real> _randCore_t;
// static _randCore_t _hw_rand_core(
//     (const float)rand_min, 
//     (const float)rand_max
// );
// Nonlinear PSO Solver
// _system_t *_hw_system_ptr = &_hw_system;
// static _randCore_t *_hw_rand_core_ptr = &_hw_rand_core;

// typedef PSO<_hw_top_real, _n_S, _maxiter, _Nh, _Nx, _n_U, _Nu> T_solver;
// T_solver my_solver(
//     _stable_zero,
//     _max_v,
//     _w0,
//     _wf,
//     _slope,
//     _c1,
//     _c2,
//     _u_min,
//     _u_max,
//     _du_max,
//     _uss
//     // ,
//     // _hw_system,
//     // _hw_rand_core_ptr
// );

/*****************************************************************************/
/**
*
* The purpose of this function is to be a wrapper for PSO as non linear solver
*
*
******************************************************************************/

int nonlinear_solver_wrapper(
    volatile float *x_curr,//[_Nx], 
    volatile float *u_curr,//[_n_U], 
    int iteration, 
    volatile float *last_best,//[_Nu*_n_U], 
    volatile float *xref,//[_Nu*_Nx], 
    // volatile float *uref,//[_n_U], 
    // volatile float *xss,//[_Nx],
    // volatile float *uss,//[_n_U], 
    
    float *new_best,//[_Nu*_n_U],
    float *J
){
#pragma HLS INTERFACE mode=s_axilite    port=iteration  bundle=control
#pragma HLS INTERFACE mode=s_axilite    port=J          bundle=control

#pragma HLS INTERFACE mode=s_axilite    port=x_curr     bundle=control
#pragma HLS INTERFACE mode=s_axilite    port=u_curr     bundle=control
#pragma HLS INTERFACE mode=s_axilite    port=last_best  bundle=control
#pragma HLS INTERFACE mode=s_axilite    port=xref       bundle=control

#pragma HLS INTERFACE mode=s_axilite    port=new_best   bundle=control

#pragma HLS INTERFACE mode=s_axilite    port=return     bundle=control

#pragma HLS INTERFACE mode=m_axi    port=x_curr     offset=slave bundle=gmem0 depth=_pragma_Nx
#pragma HLS INTERFACE mode=m_axi    port=u_curr     offset=slave bundle=gmem0 depth=_pragma_n_U
#pragma HLS INTERFACE mode=m_axi    port=last_best  offset=slave bundle=gmem0 depth=_pragma_n_U*_pragma_Nu
#pragma HLS INTERFACE mode=m_axi    port=xref       offset=slave bundle=gmem0 depth=_pragma_Nx*_pragma_Nh
#pragma HLS INTERFACE mode=m_axi    port=new_best   offset=slave bundle=gmem0 depth=_pragma_n_U*_pragma_Nu

    int iterations;
	_hw_top_real my_x_curr[_Nx] ;
	_hw_top_real my_u_curr[_n_U] ;
	_hw_top_real my_last_best[_n_U*_Nu] ;
	_hw_top_real my_xref[_Nx*_Nh];
	// _hw_top_real my_uref[_n_U] ;
	// _hw_top_real my_xss[_Nx];
	// _hw_top_real my_uss[_n_U] ;

	_hw_top_real my_new_best[_Nu*_n_U];
	_hw_top_real my_J;

//#pragma HLS bind_storage variable=local_control_guess type=FIFO impl=LUTRAM

    memcpy_loop_rolled<_hw_top_real, float, _Nx>(my_x_curr, (float *)x_curr );
    memcpy_loop_rolled<_hw_top_real, float, _n_U>(my_u_curr, (float *)u_curr );
    memcpy_loop_rolled<_hw_top_real, float, _Nu*_n_U>(my_last_best, (float *)last_best );
    memcpy_loop_rolled<_hw_top_real, float, _Nu*_Nx>(my_xref, (float *)xref );
    // memcpy_loop_rolled<_hw_top_real, float, _n_U>(my_uref, (float *)uref );
    // memcpy_loop_rolled<_hw_top_real, float, _Nx>(my_xss, (float *)xss );
    // memcpy_loop_rolled<_hw_top_real, float, _n_U>(my_uss, (float *)uss );
    
    iterations = execute(    
        (_hw_top_real *)my_x_curr, 
        (_hw_top_real *)my_u_curr, 
        iteration, 
        (_hw_top_real *)my_last_best, 
        (_hw_top_real *)my_xref,
        // my_uref, 
        // my_xss,
        (_hw_top_real *)my_new_best,
        &my_J
    );

    memcpy_loop_rolled<float, _hw_top_real, _Nu*_n_U>(new_best, (_hw_top_real *)my_new_best );
    J[0] = (float) my_J;
    return iterations;

}
