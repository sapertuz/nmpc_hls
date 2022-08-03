#ifndef HLS_PSO_HPP
#define HLS_PSO_HPP

//#define DEBUG_HLS
#include <stdint.h>

#include "config.hpp"
// #include "hls_system.hpp"
// #include "hls_pseudorand.hpp"
// #include "aux_functions.hpp"

#define H_MAX 60000.0f


// #if _Parametrization > 0
//     #define _Nn     _Parametrization
//     #define _pso_Nuu    _N
//     #define _pso_n_U     1
// #else
//     #define _Nn     _N
//     #define _pso_Nuu    _pso_Nu
//     #define _pso_n_U     _pso_n_U
// #endif 

#define _pso_hw_real 	_hw_top_real
#define _pso_n_S 		_n_S
#define _pso_maxiter 	_maxiter
#define _pso_Nh 		_Nh
#define _pso_Nx 		_Nx
#define _pso_n_U 		_n_U
#define _pso_Nu 		_Nu

// Pseudo random generator variables
// namespace pso_solver
// {

// ---------------------------------------------------
int execute(
	volatile _pso_hw_real x_curr[_pso_Nx], 
	volatile _pso_hw_real u_curr[_pso_n_U], 
	int iteration, 
	volatile _pso_hw_real last_best[_pso_n_U*_pso_Nu], 
	volatile _pso_hw_real xref[_pso_Nx*_pso_Nh],
	// volatile _pso_hw_real uref[_pso_n_U], 
	// volatile _pso_hw_real xss[_pso_Nx],
	//volatile _pso_hw_real uss[_pso_n_U], 
	
	_pso_hw_real new_best[_pso_Nu*_pso_n_U],
	_pso_hw_real J[1]
);

// Misc Functions
// ---------------------------------------------------
/*
_pso_hw_real rand_real(unsigned core){
#if (defined(__SYNTHESIS__) || defined(SYNTH_RAND))
	_pso_hw_real return_value = randGen.rand_num(core);
#else
	_pso_hw_real return_value = (_pso_hw_real)rand()/(_pso_hw_real)RAND_MAX;
#endif
	return return_value;
}
*/
// ---------------------------------------------------

void initializeParticles_set(
	// Top Level inputs (FIFO Streams)
	volatile _pso_hw_real *u_curr,
	volatile _pso_hw_real *x_curr,
	volatile _pso_hw_real *xref,
	volatile _pso_hw_real *last_best,

	// Local memories for system constraints created now
	_pso_hw_real u_curr_local[_pso_n_U],
	_pso_hw_real x_curr_local[_pso_Nx],
	_pso_hw_real xref_local[_pso_Nx*_pso_Nh],

	// Particle Data Memory
	_pso_hw_real local_x[_pso_n_S * _pso_Nu*_pso_n_U],
	_pso_hw_real local_y[_pso_n_S * _pso_Nu*_pso_n_U],
	_pso_hw_real local_v[_pso_n_S * _pso_Nu*_pso_n_U],

	// Particle Variables
	_pso_hw_real f_ind_local[_pso_n_S],

	// Local memories for system constraints created now
	_pso_hw_real local_du_min[_pso_n_U],
	_pso_hw_real local_du_max[_pso_n_U],
	_pso_hw_real local_u_min[_pso_n_U],
	_pso_hw_real local_u_max[_pso_n_U],	
	_pso_hw_real local_uss[_pso_n_U]
	
#ifdef __SYNTHESIS__
	,_rand_real_stream &__rand_port
#endif
);

// ---------------------------------------------------

void rand_real(_pso_hw_real &rand_value);

// Init Functions
// ---------------------------------------------------
// void initializeConstrains(
// 	volatile _pso_hw_real *u_curr,
// 	_pso_hw_real *local_du_max,
// 	_pso_hw_real *local_du_min,
// 	_pso_hw_real *local_u_max,
// 	_pso_hw_real *local_u_min

	// _pso_hw_real *local_x_max_first,
	// _pso_hw_real *local_x_min_first
// );
// ---------------------------------------------------
void calculate_du_min(
	_pso_hw_real *_local_du_min,
	_pso_hw_real *_local_du_max
);
// ---------------------------------------------------
void equalizeParticles(
	_pso_hw_real *local_x,
	_pso_hw_real iteration
);
// ---------------------------------------------------
_pso_hw_real verifyControlConstrains(
	_pso_hw_real local_u_max,
	_pso_hw_real local_u_min,
	_pso_hw_real value
);
#ifdef PSO_CANON
// ---------------------------------------------------
void detectInvalidParticles(
	int iter, 
	int best_pos
	
//	_pso_hw_real _x[][_pso_Nu*_pso_n_U]
	// _pso_hw_real *_x
);
#endif

// ---------------------------------------------------
void initializeParticlesWithDuConstrains(
	_pso_hw_real *u_curr,
	_pso_hw_real *local_du_max,
	
	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min,

	_pso_hw_real *local_x,
	_pso_hw_real *local_y,
	_pso_hw_real *local_v
	
	// ,bool *valid_particle

#ifdef __SYNTHESIS__
	,_rand_real_stream &__rand_port
#endif
);

#ifdef PSO_CANON
// ---------------------------------------------------
void initializeParticles(
    // _pso_hw_real _x[][_pso_Nu*_pso_n_U], 
    // _pso_hw_real _y[][_pso_Nu*_pso_n_U], 
    // _pso_hw_real _v[][_pso_Nu*_pso_n_U]
	// _pso_hw_real **_x, 
    // _pso_hw_real **_y, 
    // _pso_hw_real **_v
);
#endif
// ---------------------------------------------------
// void initializeLastBestKPSO(
// 	_pso_hw_real volatile *last_best,

// 	_pso_hw_real *local_x_max_first,
// 	_pso_hw_real *local_x_min_first,
// 	_pso_hw_real *local_x,

// 	int index
// );
// ---------------------------------------------------
void initializeStableZero(
	_pso_hw_real *uss_local, 
	_pso_hw_real *u_curr,

	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min,
	_pso_hw_real *local_du_max,
	_pso_hw_real *local_du_min,
	
	_pso_hw_real *local_x,
	
	int index
);
// Workflow Functions
// ---------------------------------------------------
void evaluateFitnessAndDetectLocalBest(
	_pso_hw_real *local_x,
	_pso_hw_real *local_y,

	_pso_hw_real *local_x_curr,//[_pso_Nx],

	_pso_hw_real *local_xref,//[_pso_Nx*_pso_Nu], 
	_pso_hw_real *local_f_ind//[_pso_n_U],
	// _pso_hw_real *local_fx//[_pso_n_U],
);
// ---------------------------------------------------
#ifdef PSO_CANON
void updateParticles(
//	_pso_hw_real global_min[],

    // _pso_hw_real _x[][_pso_Nu*_pso_n_U],
    // _pso_hw_real _y[][_pso_Nu*_pso_n_U],
    // _pso_hw_real _v[][_pso_Nu*_pso_n_U]
	// _pso_hw_real **_x,
	// _pso_hw_real **_y,
	// _pso_hw_real **_v
);
#endif
// ---------------------------------------------------
void updateParticlesWithDuConstrains(
	_pso_hw_real *local_u_curr,

	// _pso_hw_real local_w,

	_pso_hw_real *local_x,
	_pso_hw_real *local_y,
	_pso_hw_real *local_v,

	_pso_hw_real *local_global_min,

	// _pso_hw_real *local_x_max_first,
	// _pso_hw_real *local_x_min_first,

	_pso_hw_real *local_du_max,
	_pso_hw_real *local_du_min,
	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min
);

// ---------------------------------------------------
void detectGlobalMinimum(
	_pso_hw_real *local_find,//[_pso_n_S], 
	_pso_hw_real *local_bestfitness,
	_pso_hw_real *local_global_min,
	int interaction,
	
	_pso_hw_real *local_y
);
// }
#endif