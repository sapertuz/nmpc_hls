#ifndef PSO_HPP
#define PSO_HPP

#include <iostream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include "config.h"
#include "system.hpp"
#include "plot.hpp"
#include "random_number.hpp"
#include "export_data.hpp"
#include "nonlinear_solver.hpp"

class PSO : public NonlinearSolver {

private: 
	int S; 	// number of particles
	//System * controled_system;

	int kpso;
	int stable_zero;

	// PSO Parameters
	int maxiter;
	_real max_v;
	_real w0; // initial weight
	_real wf; // final weight
	_real w;
	_real slope;
	_real c1; // cognitive coefficient
	_real c2; // social coefficient

	_real * x_min;
	_real * x_max;
	_real * du_max;
	_real * x_max_first;
	_real * x_min_first;

	_real ini_v; // initial velocity
	_real threshold;
	int stop_criteria;

	// Particles
	_real ** x;
	_real ** y;
	_real ** v;
	int * valid_particle;
	int number_of_active_particles;

	_real * f_ind;
	_real * fx;
	_real * bestfitness;
	_real * global_min;
	_real * previous_control;

	// Controlled System Configuration
	int N;
	int Nu;
	int Nc;

	_real * u_from_parameters;

	RandomNumber * random;
	//ExportData * exporter;

public:
	PSO(System * controled_system, std::string config_file);
	~PSO();
	int execute(_real * u_curr, int iteration, _real * last_best, _real * xref, _real * uref, _real * xss, _real * uss, _real * J, Plot * graph);

	_real * get_x_max_first();
	_real * get_x_min_first();
	_real * get_f_ind();
	_real * get_global_min();

	int getS();
	int getKPSO();
	int getStableZero();
	int getMaxiter();
	_real getMaxV();
	int getStopCriteria();
	_real getC1();
	_real getC2();

private:
	void initializeConstrains(_real * u_curr);
	void initializeParticles();
	void initializeParticlesWithDuConstrains(_real * u_curr, _real * uss);
	void initializeLastBestKPSO(_real * last_best);
	void initializeLastBestKPSOParameters(_real * last_best);
	void initializeStableZero(_real * uss, _real * u_curr, int index);
    void initializeStableZeroParameters(_real * uss, int index);
	void initializeBestLocalFitness(void);
	void evaluateFitnessAndDetectLocalBest(_real ** x,  _real * xref, _real * uref, _real * xss, _real * uss);
	void evaluateFitnessAndDetectLocalBestParameters(_real ** x, _real * xref, _real * uref, _real * xss, _real * uss);
	int detectGlobalMinimum(int iter);
	void updateParticles(_real ** x, _real ** y, _real ** v);
	void updateParticlesWithDuConstrains(_real ** x, _real ** y, _real ** v, int test);
	void equalizeParticlesVant();
	void verifyControlConstrains(_real * value, int pos);
	void detectInvalidParticles(int iter, int best_pos);

};

#endif
