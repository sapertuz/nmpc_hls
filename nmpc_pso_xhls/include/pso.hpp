#ifndef PSO_HPP
#define PSO_HPP

#define DEBUG_HLS

#include <stdio.h>
#include <stdlib.h>
#include "config.hpp"
#include "system.hpp"
#include "nonlinear_solver.hpp"

#if _Parametrization > 0
    #define _Nn     _Parametrization
    #define _Nuu    _N
    #define _Nc     1
#else
    #define _Nn     _N
    #define _Nuu    _Nu
    #define _Nc     _n_U
#endif 


class PSO : public NonlinearSolver {

private: 
	const _uchar S = _S; 	// number of particles
	//System * controled_system;

	const _uchar kpso = _KPSO;
	const _uchar stable_zero = _stable_zero;

	// PSO Parameters
	const _uchar maxiter = _maxiter;
	_real max_v;
	_real w0; // initial weight
	_real wf; // final weight
	_real w;
	_real slope;
	_real c1; // cognitive coefficient
	_real c2; // social coefficient
	
	_real threshold;
	int stop_criteria;

	_real x_min[_Nn];
	_real x_max[_Nn];
	_real du_max[_n_U];
	_real x_max_first[_n_U];
	_real x_min_first[_n_U];

	_real ini_v; // initial velocity

	// Particles
	_real x[_S][_Nu*_n_U];
	_real y[_S][_Nu*_n_U];
	_real v[_S][_Nu*_n_U];
	int valid_particle[_S];
	int number_of_active_particles;

	_real f_ind[_S];
	_real fx[_S];
	_real bestfitness[_maxiter];
	_real global_min[_Nu*_n_U];
	_real previous_control[_Parametrization];

	// Controlled System Configuration
	int N;
	int Nu;
	int Nc;

	_real u_from_parameters[_N*_n_U];

	//ExportData * exporter;

public:
	PSO();
	~PSO();
	int execute(
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
	);

	int getS();
	int getKPSO();
	int getStableZero();
	int getMaxiter();
	_real getMaxV();
	int getStopCriteria();
	_real getC1();
	_real getC2();

private:
	_real rand_real();
	_real min_array(
		_real array[_S], 
		int * pos
	);
	void initializeConstrains(
		_real u_curr[_n_U]
	);
	void initializeBestLocalFitness(void);
	int detectGlobalMinimum(
		int iter
	);
	void equalizeParticlesVant();
	void verifyControlConstrains(
		_real value[1], 
		int pos
	);
	void detectInvalidParticles(
		int iter, 
		int best_pos
	);

	void initializeParticlesWithDuConstrains(
		_real u_curr[_n_U], 
		_real uss[_n_U]
	);
	void initializeLastBestKPSO(
		_real last_best[_Nu*_n_U]
	);
	void initializeStableZero(
		_real uss[_n_U], 
		_real u_curr[_n_U], 
		int index
	);
    void evaluateFitnessAndDetectLocalBest(
		_real x[_S][_n_U*_Nu],
		_real xref[_Nx*_Nu], 
		_real uref[_n_U],
		_real xss[_Nx], 
		_real uss[_n_U]
	);
	void updateParticlesWithDuConstrains(
		_real x[_S][_Nu*_n_U], 
		_real y[_S][_Nu*_n_U], 
		_real v[_S][_Nu*_n_U]
	);

	void initializeParticles();	
	void updateParticles(
		_real x[_S][_Nu*_n_U], 
		_real y[_S][_Nu*_n_U], 
		_real v[_S][_Nu*_n_U]
	);
	
};

#endif
