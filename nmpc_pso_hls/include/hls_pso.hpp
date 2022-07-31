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
namespace pso_solver
{
const unsigned int n_S = _pso_n_S; 	// number of particles
//System * controled_system;

// const int kpso = _KPSO;
// const int stable_zero;
const uint8_t particle_last_best = 0;
const uint8_t particle_stable_zero = 1;

// PSO Parameters
// const unsigned int	_pso_maxiter = _pso_maxiter;
const _pso_hw_real 	_pso_max_v = _max_v;
const _pso_hw_real 	_pso_min_v = -_max_v;
const _pso_hw_real 	_pso_w0 = _w0; // initial weight
const _pso_hw_real	_pso_wf = _wf; // final weight
const _pso_hw_real	_pso_slope = _slope;
const _pso_hw_real	_pso_c1 = _c1; // cognitive coefficient
const _pso_hw_real	_pso_c2 = _c2; // social coefficient


const _pso_hw_real _pso_init_v = _max_v*0.1; // initial velocity

// Controlled System Configuration
const unsigned int Nh = _pso_Nh;
const unsigned int Nx = _pso_Nx;
const unsigned int n_U = _pso_n_U;
const unsigned int Nu = _pso_Nu;
const unsigned int part_S = Nu*n_U;
// _pso_hw_real u_from_parameters[_N*_pso_n_U];

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


void rand_real(_pso_hw_real &rand_value);

// ---------------------------------------------------
_pso_hw_real min_array(
	_pso_hw_real *array,//[_pso_n_S], 
	int * pos
);

// Init Functions
// ---------------------------------------------------
void initializeConstrains(
	volatile _pso_hw_real *u_curr
);
// ---------------------------------------------------
void calculate_du_min(
);
// ---------------------------------------------------
void initializeBestLocalFitness(void);
// ---------------------------------------------------
int detectGlobalMinimum(
	int iter
);
// ---------------------------------------------------
void equalizeParticles(
//	_pso_hw_real _x[][_pso_Nu*_pso_n_U],
	// _pso_hw_real **_x,
	_pso_hw_real iteration
);
// ---------------------------------------------------
_pso_hw_real verifyControlConstrains(
	_pso_hw_real value, 
	int pos
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
// Workflow Functions
// ---------------------------------------------------
void initializeParticlesWithDuConstrains(
	volatile _pso_hw_real *u_curr 
	//volatile _pso_hw_real *uss,

	// _pso_hw_real _x[][_pso_Nu*_pso_n_U], 
	// _pso_hw_real _y[][_pso_Nu*_pso_n_U], 
	// _pso_hw_real _v[][_pso_Nu*_pso_n_U]
	// _pso_hw_real ** _x,
	// _pso_hw_real ** _y,
	// _pso_hw_real ** _v
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
void initializeLastBestKPSO(
	_pso_hw_real volatile *last_best,
	// _pso_hw_real _x[][_pso_Nu*_pso_n_U]
	// _pso_hw_real **_x,
	uint8_t index
);
// ---------------------------------------------------
void initializeStableZero(
//	_pso_hw_real *uss, 
	_pso_hw_real *u_curr, 	
	// _pso_hw_real _x[][_pso_Nu*_pso_n_U]
	// _pso_hw_real **_x,
	uint8_t index
);
// ---------------------------------------------------
void evaluateFitnessAndDetectLocalBest(
	// _pso_hw_real _x[][_pso_Nu*_pso_n_U], 
	// _pso_hw_real _y[][_pso_Nu*_pso_n_U],
	// _pso_hw_real **_x,
	// _pso_hw_real **_y,

	_pso_hw_real *x_curr_local,//[_pso_Nx],

	_pso_hw_real *xref//[_pso_Nx*_pso_Nu], 
	// _pso_hw_real *uref,//[_pso_n_U],
	// _pso_hw_real *xss//[_pso_Nx]
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
	_pso_hw_real *u_curr

//	_pso_hw_real global_min[],

	// _pso_hw_real *_x[][_pso_Nu*_pso_n_U], 
	// _pso_hw_real *_y[][_pso_Nu*_pso_n_U], 
	// _pso_hw_real *_v[][_pso_Nu*_pso_n_U]
	// _pso_hw_real **_x,
	// _pso_hw_real **_y,
	// _pso_hw_real **_v
);

}
#endif