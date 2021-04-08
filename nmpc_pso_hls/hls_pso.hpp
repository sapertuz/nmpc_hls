#ifndef HLS_PSO_HPP
#define HLS_PSO_HPP

#define DEBUG_HLS

#include "config.hpp"
//#include "system.hpp"

#if _Parametrization > 0
    #define _Nn     _Parametrization
    #define _Nuu    _N
    #define _Nc     1
#else
    #define _Nn     _N
    #define _Nuu    _Nu
    #define _Nc     _n_U
#endif 

// Pseudo random generator variables
const int a = (1103515245);
const int c = (12345);
const int m = (1<<31);
short rnd_seed = 9727;

// DEBUG

struct HLS_PSO {
	const int S = _S; 	// number of particles
	//System * controled_system;

	const int kpso = _KPSO;
	const int stable_zero = _stable_zero;

	// PSO Parameters
	const int 		maxiter = _maxiter;
	_hw_real 	max_v = _max_v;
	_hw_real 	w0 = _w0; // initial weight
	_hw_real	wf = _wf; // final weight
	_hw_real	slope = (wf-w0)/maxiter;
	_hw_real	c1 = _c1; // cognitive coefficient
	_hw_real	c2 = _c2; // social coefficient
	
	_hw_real w = w0;
	const _hw_real slope_init = (wf-w0)/maxiter;
	
	_hw_real threshold;
	int stop_criteria;

	_hw_real x_min[_Nn];
	_hw_real x_max[_Nn];
	_hw_real du_max[_n_U];
	_hw_real x_max_first[_n_U];
	_hw_real x_min_first[_n_U];

	const _hw_real ini_v = max_v*0.01; // initial velocity

	// Particles
	// _hw_real x[_S][_Nu*_n_U];
	// _hw_real y[_S][_Nu*_n_U];
	// _hw_real v[_S][_Nu*_n_U];
	int valid_particle[_S];
	int number_of_active_particles;

	_hw_real f_ind[_S];
	_hw_real fx[_S];
	// _hw_real bestfitness[_maxiter];
	// _hw_real global_min[_Nu*_n_U];
	_hw_real previous_control[_Parametrization];

	// Controlled System Configuration
	int N;
	int Nu;
	int Nc;

	// _hw_real u_from_parameters[_N*_n_U];

	//ExportData * exporter;
};

// Init Execute
void hls_pso__init(HLS_PSO * myself);

int hls_pso__execute(
	HLS_PSO * myself,
	_hw_real x_curr[_Nx], 
	_hw_real u_curr[_n_U], 
	int iteration, 
	_hw_real last_best[_Nu*_n_U], 
	_hw_real xref[_Nu*_Nx],
	_hw_real uref[_n_U], 
	_hw_real xss[_Nx],
	_hw_real uss[_n_U], 
	_hw_real new_best[_Nu*_n_U],
	_hw_real * J
);

// Misc Functions
_hw_real hls_pso__rand_real();
_hw_real hls_pso__min_array(
	_hw_real array[_S], 
	int * pos
);
// Init Functions
void hls_pso__initializeConstrains(
	_hw_real u_curr[_n_U],
    _hw_real x_max_first[],
    _hw_real x_min_first[]
);
void hls_pso__initializeBestLocalFitness(void);
int hls_pso__detectGlobalMinimum(
	int iter
);
void hls_pso__equalizeParticles(
	_hw_real _x[_S][_Nu*_n_U]
);
void hls_pso__verifyControlConstrains(
	_hw_real value[1], 
	int pos
);
void hls_pso__detectInvalidParticles(
	int iter, 
	int best_pos,
	
	_hw_real _x[_S][_Nu*_n_U]
);

// Workflow Functions
void hls_pso__initializeParticlesWithDuConstrains(
	_hw_real u_curr[_n_U], 
	_hw_real uss[_n_U],

	_hw_real _x[_S][_Nu*_n_U], 
	_hw_real _y[_S][_Nu*_n_U], 
	_hw_real _v[_S][_Nu*_n_U]
);
void hls_pso__initializeLastBestKPSO(
	_hw_real last_best[_Nu*_n_U],
	
	_hw_real _x[_S][_Nu*_n_U]
);
void hls_pso__initializeStableZero(
	_hw_real uss[_n_U], 
	_hw_real u_curr[_n_U], 
	int index,
	
	_hw_real _x[_S][_Nu*_n_U]
);
void hls_pso__evaluateFitnessAndDetectLocalBest(
	_hw_real _x[_S][_Nu*_n_U], 
	_hw_real _y[_S][_Nu*_n_U],

	_hw_real xref[_Nx*_Nu], 
	_hw_real uref[_n_U],
	_hw_real xss[_Nx], 
	_hw_real uss[_n_U]
);
void hls_pso__updateParticlesWithDuConstrains(
	_hw_real _x[_S][_Nu*_n_U], 
	_hw_real _y[_S][_Nu*_n_U], 
	_hw_real _v[_S][_Nu*_n_U]
);

#endif