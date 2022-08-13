#ifndef HLS_PSO_HPP
#define HLS_PSO_HPP

//#define DEBUG_HLS
#include <stdint.h>

#include "config.hpp"
// #include "hls_system.hpp"
// #include "hls_pseudorand.hpp"
// #include "aux_functions.hpp"

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
	_pso_hw_real x_curr[_pso_Nx], 
	_pso_hw_real u_curr[_pso_n_U], 
	int iteration, 
	_pso_hw_real last_best[_pso_n_U*_pso_Nu], 
	_pso_hw_real xref[_pso_Nx*_pso_Nh],
	// volatile _pso_hw_real uref[_pso_n_U], 
	// volatile _pso_hw_real xss[_pso_Nx],
	//volatile _pso_hw_real uss[_pso_n_U], 
	
	_pso_hw_real new_best[_pso_Nu*_pso_n_U],
	_pso_hw_real J[1]
);

// Top Functions

// ---------------------------------------------------

void initializeParticles_set(
	// Top Level inputs (FIFO Streams)
#ifndef __SYNTHESIS__
	_pso_hw_real u_curr[_pso_n_U],
	_pso_hw_real x_curr[_pso_Nx],
	_pso_hw_real xref[_pso_Nx*_pso_Nh],
#endif
	_pso_hw_real last_best[_pso_Nu*_pso_n_U],

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

// Aux Functions

// ---------------------------------------------------

void rand_real(_pso_hw_real &rand_value);

// Init Functions

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
// ---------------------------------------------------
void detectInvalidParticles(
	int iter, 
	int best_pos
	
);

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

// ---------------------------------------------------
void initializeParticles(
	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min,

	_pso_hw_real *local_x,
	_pso_hw_real *local_y,
	_pso_hw_real *local_v
);

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

void updateParticles(
	_pso_hw_real *local_x,
	_pso_hw_real *local_y,
	_pso_hw_real *local_v,

	_pso_hw_real *local_global_min,

	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min
);

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

#ifdef __SYNTHESIS__
	,_rand_real_stream &__rand_port
#endif
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