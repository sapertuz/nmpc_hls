#include "hls_pso.hpp"

#include "hls_system.hpp"
#include "hls_pseudorand.hpp"
#include "aux_functions.hpp"

using namespace pso_solver;

const unsigned int n_S = _pso_n_S; 	// number of particles
//System * controled_system;

// const int kpso = _KPSO;
// const int stable_zero;
const int particle_last_best = 0;
const int particle_stable_zero = 1;

// PSO Parameters
// const unsigned int	_pso_maxiter = _pso_maxiter;
const _pso_hw_real 	_pso_max_v = _max_v;
const _pso_hw_real 	_pso_min_v = -_max_v;
const _pso_hw_real 	_pso_w0 = _w0; // initial weight
const _pso_hw_real	_pso_wf = _wf; // final weight
const _pso_hw_real	_pso_slope = _slope;
const _pso_hw_real	_pso_c1 = _c1; // cognitive coefficient
const _pso_hw_real	_pso_c2 = _c2; // social coefficient

const _pso_hw_real *_pso_u_max(_u_max);
const _pso_hw_real *_pso_u_min(_u_min);
const _pso_hw_real *_pso_uss_local(_uss);//[_pso_n_U];

// const _pso_hw_real *_pso_u_min(_u_min);//[_pso_n_U];
// const _pso_hw_real *_pso_u_max(_u_max);//[_pso_n_U];

const _pso_hw_real _pso_init_v = _max_v*0.1; // initial velocity

// Controlled System Configuration
const unsigned int Nh = _pso_Nh;
const unsigned int Nx = _pso_Nx;
const unsigned int n_U = _pso_n_U;
const unsigned int Nu = _pso_Nu;
const unsigned int part_S = Nu*n_U;
// _pso_hw_real u_from_parameters[_N*_pso_n_U];


// Pseudo Random Core Generator
typedef pseudoRand_gen<_pso_hw_real> _pso_randCore_t;
// _pso_randCore_t randGen(
// 	(const float)pso_rand_min, 
// 	(const float)pso_rand_max
// );
_pso_randCore_t rand_core;

// ---------------------------------------------------

int pso_solver::execute(
	volatile _pso_hw_real *x_curr,//[_pso_Nx], 
	volatile _pso_hw_real *u_curr,//[_pso_n_U], 
	int iteration, 
	volatile _pso_hw_real *last_best,//[_pso_n_U*_pso_Nu], 
	volatile _pso_hw_real *xref,//[_pso_Nx*_pso_Nh],
	// volatile _pso_hw_real uref[_pso_n_U], 
	// volatile _pso_hw_real xss[_pso_Nx],
	//volatile _pso_hw_real uss[_pso_n_U], 
	
	_pso_hw_real *new_best,//[_pso_Nu*_pso_n_U],
	_pso_hw_real *J//[1]
){

	_pso_hw_real w;

	_pso_hw_real threshold;
	int stop_criteria;

	// _pso_hw_real x_max_first[_pso_n_U];
	// _pso_hw_real x_min_first[_pso_n_U];


	_pso_hw_real *_pso_du_max((_pso_hw_real *)_du_max);//[_pso_n_U];
	_pso_hw_real _pso_du_min[_pso_n_U];//[_pso_n_U];

// Particles
#ifdef __SYNTHESIS__
	_pso_hw_real x[_pso_n_S * _pso_Nu*_pso_n_U];
	_pso_hw_real y[_pso_n_S * _pso_Nu*_pso_n_U];
	_pso_hw_real v[_pso_n_S * _pso_Nu*_pso_n_U];
#else
	//_pso_hw_real bestfitness[_pso_maxiter];
	//_pso_hw_real global_min[Nu*n_U];
	_pso_hw_real *x;
	_pso_hw_real *y;
	_pso_hw_real *v;
#endif
	bool valid_particle[_pso_n_S];
	int number_of_active_particles;

	_pso_hw_real f_ind[_pso_n_S];
	_pso_hw_real fx[_pso_n_S];
	_pso_hw_real bestfitness[_pso_maxiter];
	_pso_hw_real global_min[_pso_Nu*_pso_n_U];
	_pso_hw_real _pso_x_min_first[_pso_n_U], _pso_x_max_first[_pso_n_U];
	//_pso_hw_real previous_control[_Parametrization];


	// _pso_system_t current_system;
	// _pso_randCore_t randGen;

	// const float pso_rand_min = 0.0f;
	// const float pso_rand_max = 1.0f;

#pragma HLS ALLOCATION operation instances=hadd limit=2
#pragma HLS ALLOCATION operation instances=hsub limit=2
#pragma HLS ALLOCATION operation instances=hmul limit=2

	// Update sensor read
    //static System * controled_system = new ModelState();
    //controled_system->setState(x_curr);

	// Particles Variables
#ifndef __SYNTHESIS__
	x = (_pso_hw_real *) malloc(n_S * Nu*n_U*sizeof(_pso_hw_real)); // alloc_matrix(n_S, Nu*n_U)) == NULL) {return -1;}
	y = (_pso_hw_real *) malloc(n_S * Nu*n_U*sizeof(_pso_hw_real)); // alloc_matrix(n_S, Nu*n_U)) == NULL) {return -1;}
	v = (_pso_hw_real *) malloc(n_S * Nu*n_U*sizeof(_pso_hw_real)); // alloc_matrix(n_S, Nu*n_U)) == NULL) {return -1;}
#endif

    number_of_active_particles = n_S;    
	int best_pos;

#ifdef __SYNTHESIS__
	_pso_hw_real u_curr_local[_pso_n_U];
#pragma HLS bind_storage variable=u_curr_local type=RAM_2P impl=LUTRAM

	_pso_hw_real x_curr_local[_pso_Nx];
#pragma HLS bind_storage variable=x_curr_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=x_curr_local type=shared

	_pso_hw_real xref_local[_pso_Nx*_pso_Nh];
#pragma HLS bind_storage variable=xref_local type=RAM_2P impl=BRAM
#pragma HLS stream variable=xref_local type=shared

#else
	_pso_hw_real *u_curr_local;
	_pso_hw_real *x_curr_local;
	_pso_hw_real *xref_local;
	// _pso_hw_real *uref_local;
	// _pso_hw_real *xss_local;
	x_curr_local = (_pso_hw_real *) malloc(Nx*sizeof(_pso_hw_real));
	u_curr_local = (_pso_hw_real *) malloc(n_U*sizeof(_pso_hw_real));
	xref_local = (_pso_hw_real *) malloc(Nx*Nh*sizeof(_pso_hw_real));
	// uref_local = (_pso_hw_real *) malloc(n_U*sizeof(_pso_hw_real));
	// xss_local = (_pso_hw_real *) malloc(Nx*sizeof(_pso_hw_real));
#endif


	pso_solver::initializeParticles_set(
		u_curr,
		x_curr,
		xref,
		last_best,

		x,y,v,

		xref_local,
		f_ind,
		w,

		u_curr_local,
		x_curr_local,

		_pso_du_min,
		_pso_du_max,
		(_pso_hw_real *)_pso_u_min,
		(_pso_hw_real *)_pso_u_max,
		_pso_x_min_first,
		_pso_x_max_first,
		(_pso_hw_real *)_pso_uss_local,
		valid_particle
	);

	// ITERATIVE PROCESS
	unsigned int k = 0;  // index of iteration

#ifdef DEBUG_PSO
	std::cout << "iter   \t";
	for (unsigned int i = 0; i < n_U; i++)
		std::cout << "u[" << i << "] \t";
	std::cout << "|       ";
	for (unsigned int i = 0; i < n_S; i++)
		std::cout << "J[" << i << "]       \t";		
	std::cout << "|";
	std::cout << std::endl;
#endif

	while (k < _pso_maxiter) {
        pso_solver::evaluateFitnessAndDetectLocalBest(
			x, 
			y, 
			
			x_curr_local, 

			xref_local,
			f_ind,
			fx
		);

        // Global Minimum detection
		// best_pos = pso_solver::detectGlobalMinimum(k);
		pso_solver::detectGlobalMinimum(
			f_ind,
			bestfitness,
			global_min,
			k,
			y
		);
		// memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nu*_pso_n_U>(global_min, (_pso_hw_real *)&y[best_pos * part_S]);

        pso_solver::updateParticlesWithDuConstrains(
			u_curr_local,
			w,

			x,
			y,
			v,

			global_min,

			_pso_x_max_first,
			_pso_x_min_first,

			_pso_du_max,
			_pso_du_min,
			(_pso_hw_real *)_pso_u_max,
			(_pso_hw_real *)_pso_u_min
		);    
#ifdef DEBUG_PSO
		std::cout << "[" << k << "] ";
		std::cout << "\t"; print_formatted_float_array(&global_min[0], n_U, 2, 6);
		std::cout << "|";
		std::cout << "\t"; print_formatted_float_array(f_ind, n_S, 2, 6); 
		std::cout << "|";
		std::cout << std::endl;
#endif
        //detectInvalidParticles(k, best_pos);

		k++;
		w+= _pso_slope;
	}


	// Return best value
    for (unsigned int i = 0; i < Nu*n_U; ++i) {
        new_best[i] = global_min[i];
    }
    
    J[0] = bestfitness[_pso_maxiter-1];

#ifndef __SYNTHESIS__
    if(x != NULL) free(x);
    if(y != NULL) free(y);
    if(v != NULL) free(v);
	if(u_curr_local != NULL) free(u_curr_local);
	if(x_curr_local != NULL) free(x_curr_local);
	if(xref_local != NULL) free(xref_local);
	// if(uref_local != NULL) free(uref_local);
	// if(xss_local != NULL) free(xss_local);
#endif
	return k;
}

void pso_solver::initializeParticles_set(
	// Top Level inputs (FIFO Streams)
	volatile _pso_hw_real *u_curr,
	volatile _pso_hw_real *x_curr,
	volatile _pso_hw_real *xref,
	volatile _pso_hw_real *last_best,

	// Particle Data Memory
	_pso_hw_real *local_x,
	_pso_hw_real *local_y,
	_pso_hw_real *local_v,

	// Particle Variables
	_pso_hw_real *xref_local,
	_pso_hw_real *f_ind_local,
	_pso_hw_real &w,

	// Local memories for system constraints created now
	_pso_hw_real *u_curr_local,
	_pso_hw_real *x_curr_local,

	// Local memories for system constraints created now
	_pso_hw_real *local_du_min,
	_pso_hw_real *local_du_max,
	_pso_hw_real *local_u_min,
	_pso_hw_real *local_u_max,
	_pso_hw_real *local_x_min_first,
	_pso_hw_real *local_x_max_first,
	_pso_hw_real *local_uss,
	
	bool *valid_particle
){

	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(u_curr_local, 	(_pso_hw_real *)u_curr);
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nx>(x_curr_local, 	(_pso_hw_real *)x_curr);
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nx*_pso_Nh>(xref_local, (_pso_hw_real *)xref);
	// memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(uref_local, 	(_pso_hw_real *)uref);
	// memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nx>(xss_local, 		(_pso_hw_real *)xss);

	w = _pso_w0;
	pso_solver::calculate_du_min(
		local_du_min,
		local_du_max
	);

	pso_solver::initializeConstrains(
		u_curr_local,
		local_du_max,
		local_du_min,
		local_u_max,
		local_u_min,

		local_x_max_first,
		local_x_min_first
	);
    pso_solver::initializeParticlesWithDuConstrains(
		u_curr_local,
		local_du_max,
		
		local_u_max,
		local_u_min,

		local_x,
		local_y,
		local_v,
		
		valid_particle
	);

	// equalizeParticles(iteration);

    pso_solver::initializeLastBestKPSO(
		last_best, 
	
		local_x_max_first,
		local_x_min_first,
		
		local_x,

		particle_last_best
	);

	pso_solver::initializeStableZero(
		local_uss,
		u_curr_local, 

		local_u_max,
		local_u_min,
		local_du_max,
		local_du_min,
		
		local_x,

		particle_stable_zero
	);

    // pso_solver::initializeBestLocalFitness();
	memset_loop<_pso_hw_real>(f_ind_local, (const _pso_hw_real)H_MAX, n_S);
}
// ---------------------------------------------------

void pso_solver::rand_real(_pso_hw_real &rand_value){
// #pragma HLS inline off
	// static _pso_randCore_t rand_core;
	_pso_hw_real rand_num_out;
#if (defined(__SYNTHESIS__) || defined(SYNTH_RAND))
	rand_core.rand_num(rand_num_out);
#else
	rand_num_out = (_pso_hw_real)rand()/(_pso_hw_real)RAND_MAX;
#endif
	rand_value = rand_num_out;
}


// ---------------------------------------------------

void  pso_solver::detectGlobalMinimum(
	_pso_hw_real *local_find,//[_pso_n_S], 
	_pso_hw_real *local_bestfitness,
	_pso_hw_real *local_global_min,
	int interaction,
	
	_pso_hw_real *local_y
){
	//[bestfitness(k), p] = min(f_ind);
#pragma HLS inline
	_pso_hw_real min = local_find[0];
	int best_pos = 0;
	for (unsigned int i = 1; i < n_S; ++i) {
        if(local_find[i] < min) {
			min = local_find[i];
			best_pos = i;
        }
	}
	local_bestfitness[interaction] = min;
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nu*_pso_n_U>(local_global_min, (_pso_hw_real *)&local_y[best_pos * part_S]);
}

// ---------------------------------------------------

void pso_solver::initializeConstrains(
	volatile _pso_hw_real *u_curr,
	_pso_hw_real *local_du_max,
	_pso_hw_real *local_du_min,
	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min,

	_pso_hw_real *local_x_max_first,
	_pso_hw_real *local_x_min_first
){
	// Initialize constrains based on current state
#pragma HLS inline
	initializeConstrains_loop: for (unsigned int i = 0; i < n_U; ++i){
		_pso_hw_real u_curr_tmp = u_curr[i];

        _pso_hw_real x_max_first_temp = u_curr_tmp + local_du_max[i]; //controled_system->acc_max[i];
		_pso_hw_real x_min_first_temp = u_curr_tmp + local_du_min[i];//controled_system->acc_min[i];        

		local_x_max_first[i] = (x_max_first_temp > local_u_max[i]) ? local_u_max[i]: x_max_first_temp;
		local_x_min_first[i] = (x_min_first_temp < local_u_min[i])? local_u_min[i]: x_min_first_temp;	
	}
}

// ---------------------------------------------------

void pso_solver::calculate_du_min(
	_pso_hw_real *_local_du_min,
	_pso_hw_real *_local_du_max
){
#pragma HLS inline
	calculate_du_min_loop: for (unsigned i = 0; i < n_U; i++)
	{
		_local_du_min[i] = -_local_du_max[i];
	}
	
}

// ---------------------------------------------------

void pso_solver::equalizeParticles(
	_pso_hw_real *local_x,
	_pso_hw_real iteration
){
	_pso_hw_real x_tmp;
	if (iteration == 0){	
    equalizeParticles_loop_S: for (unsigned int i = 0; i < n_S; ++i) {
		equalizeParticles_loop_N:for (unsigned int k = 1; k < Nu; ++k) {
	 	   	int idx1 = k*n_U;
	 	   	int idx2 = 0;
			equalizeParticles_loop_x:for (unsigned int j = 0; j < n_U; ++j) {
				x_tmp = local_x[i*part_S + idx2];
                local_x[i*part_S + idx1] = x_tmp;
				idx1++;
				idx2++;
            }
        }
    }
	}
}

// ---------------------------------------------------

_pso_hw_real pso_solver::verifyControlConstrains(
	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min,
	_pso_hw_real value, 
	int pos
){
#pragma HLS inline
	_pso_hw_real return_value = value;

	return_value = (return_value > local_u_max[pos]) ? local_u_max[pos] : return_value;
	return_value = (return_value < local_u_min[pos]) ? local_u_min[pos] : return_value;

	return return_value;
}

// ---------------------------------------------------
#ifdef PSO_CANON

void pso_solver::detectInvalidParticles(
	int iter, 
	int best_pos
	
//	_pso_hw_real _x[][_pso_Nu*_pso_n_U]
	// _pso_hw_real *_x
){
	for (unsigned int i = 0; i < n_S; ++i) {
        if(n_S != best_pos) {
            if((std::abs(fx[i]-bestfitness[iter])/bestfitness[iter]) < 0.0001){
//                int count = 0;  
                int count = Nh*n_U;
                for (unsigned int j = 0; j < Nh*n_U; ++j) {
                    if((std::abs(x[i*part_S + j]-global_min[j])/global_min[j]) < 0.0001){
                        count--;
                    }
                }
                if((valid_particle[i] == 1) && (count == 0)){
                    valid_particle[i] = 0;
                    number_of_active_particles--;
                }
            }
        }
    }

}
#endif


// ---------------------------------------------------

void pso_solver::initializeParticlesWithDuConstrains(
	_pso_hw_real *u_curr,

	_pso_hw_real *local_du_max,
	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min,
	
	_pso_hw_real *local_x,
	_pso_hw_real *local_y,
	_pso_hw_real *local_v,
	
	bool *valid_particle
){
// INITIALIZATION OF PARTICLES WITH Delta u CONTRAINS FOR MPC
//#pragma HLS inline
#pragma HLS ALLOCATION instances=rand_real limit=1 function
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation

	_pso_hw_real x_ant[_pso_n_U];
	_pso_hw_real x_ant_tmp[_pso_n_U];

	initializeS_du_loop_S:for (unsigned int i = 0; i < n_S; ++i) {
		memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant_tmp, (_pso_hw_real *)u_curr);
		initializeS_du_loop_N:for (unsigned int j = 0; j < Nu; ++j) {
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant, (_pso_hw_real *)x_ant_tmp);
			initializeS_du_loop_x:for (unsigned int k = 0; k < n_U; ++k) {
// #pragma HLS pipeline II=11 rewind
				int idx = j*n_U + k;
				_pso_hw_real rand_tmp;
				_pso_hw_real local_du_max_k = local_du_max[k];
				pso_solver::rand_real(rand_tmp);
				_pso_hw_real du_tmp = rand_tmp*(_pso_hw_real)2.0*local_du_max_k - local_du_max_k;
		        _pso_hw_real x_tmp = x_ant[k] + du_tmp; //random->read(); 
		        // _pso_hw_real x_tmp = ( x_ant[k] + (-du_max[k]) ) + ((_pso_hw_real)2.0*du_max[k]) * rand_tmp; //random->read(); 
				x_tmp = pso_solver::verifyControlConstrains(
					local_u_max,
					local_u_min,
					x_tmp, 
					k
				);
                local_x[i*part_S + idx] = x_tmp;
				local_y[i*part_S + idx] = x_tmp; // uss_local[k];
		        local_v[i*part_S + idx] = _pso_init_v;
                valid_particle[i] = true;
				x_ant_tmp[k] = x_tmp;
				idx++;
		    }
		}
	}
}

#ifdef PSO_CANON

// ---------------------------------------------------

void  pso_solver::initializeParticles(
    // _pso_hw_real _x[][_pso_Nu*_pso_n_U], 
    // _pso_hw_real _y[][_pso_Nu*_pso_n_U], 
    // _pso_hw_real _v[][_pso_Nu*_pso_n_U]
	// _pso_hw_real **_x, 
    // _pso_hw_real **_y, 
    // _pso_hw_real **_v
) {
// GENERAL PSO INITIALIZATION OF PARTICLES
	_pso_hw_real d_max = 0.2;
    _pso_hw_real d2_max = 0.4;
	_pso_hw_real x_ant[n_U];
	for (unsigned int i = 0; i < n_S; ++i) {
		for (unsigned int j = 0; j < Nu; ++j) {
			for (unsigned int k = 0; k < n_U; ++k) {
				int idx = j*n_U + k;
				_pso_hw_real rand_tmp; 
				rand_real(rand_tmp);
				_pso_hw_real x_tmp = (j == 0)? 
					u_min[k] + (u_max[k]-u_min[k]) * rand_tmp: //random->read();
					x_ant[k] - d_max + (d2_max) * rand_tmp; //random->read();
				x[i*part_S + idx] = x_tmp; 
				y[i*part_S + idx] = H_MAX;
				v[i*part_S + idx] = init_v;  
				x_ant[k] = x_tmp;
			}
        }
        valid_particle[i] = 1;
    }
}
#endif

// ---------------------------------------------------

void pso_solver::initializeLastBestKPSO(
	_pso_hw_real volatile *last_best,

	_pso_hw_real *local_x_max_first,
	_pso_hw_real *local_x_min_first,
	_pso_hw_real *local_x,

	int index
){
// #pragma HLS inline
	// Uses KPSO (best last position)

	for (unsigned int i = 0; i < n_U; i++){
		_pso_hw_real last_best_tmp = last_best[i];
		_pso_hw_real x_i;
		x_i = (last_best_tmp > local_x_max_first[i]) ? local_x_max_first[i] : last_best_tmp;
		x_i = (last_best_tmp < local_x_min_first[i]) ? local_x_max_first[i] : last_best_tmp;
		local_x[index*part_S + i] = x_i;
	}
	
	for (unsigned int i = n_U; i < (Nu*n_U); i++) {
		local_x[index*part_S + i] = last_best[i];			
	}
}

// ---------------------------------------------------

void pso_solver::initializeStableZero(
	_pso_hw_real *uss_local, 
	_pso_hw_real *u_curr,

	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min,
	_pso_hw_real *local_du_max,
	_pso_hw_real *local_du_min,

	_pso_hw_real *local_x,
	
	int index
){
// #pragma HLS inline
// Starts one particle with all Zeros for 
// stable response after equilibrium is reached
    int idx;
	_pso_hw_real x_ant[_pso_n_U];
	_pso_hw_real x_ant_tmp[_pso_n_U];
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant_tmp, (_pso_hw_real*)u_curr);
   	for (unsigned int k = 0; k < Nu; ++k){
		memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant, (_pso_hw_real*)x_ant_tmp);
        for (unsigned int i = 0; i < n_U; ++i) {
        	idx = k*n_U + i;
		    _pso_hw_real x_tmp = uss_local[i];
			// _pso_hw_real comp_tmp = (k == 0) ? x_tmp - u_curr[i] : x_tmp - x_ant[i];
			_pso_hw_real comp_tmp = x_tmp - x_ant[i] ;
			
			x_tmp = (comp_tmp > local_du_max[i]) ? local_du_max[i] : x_tmp;
			x_tmp = (comp_tmp < local_du_min[i]) ? local_du_min[i] : x_tmp;

            local_x[index*part_S + idx] = pso_solver::verifyControlConstrains(
				local_u_max,
				local_u_min,
				x_tmp, 
				i
			);
			x_ant_tmp[i] = x_tmp;
			idx++;
        }
   	}
}

// ---------------------------------------------------

void pso_solver::evaluateFitnessAndDetectLocalBest(
	_pso_hw_real *local_x,
	_pso_hw_real *local_y,

	_pso_hw_real *local_x_curr,//[_pso_Nx],

	_pso_hw_real *local_xref,//[_pso_Nx*_pso_Nu], 
	_pso_hw_real *local_f_ind,//[_pso_n_U],
	_pso_hw_real *local_fx//[_pso_n_U],
){

// #pragma HLS interface ap_bus  port=_x
// #pragma HLS interface ap_bus  port=xref
// #pragma HLS interface ap_fifo port=uref
// #pragma HLS interface ap_fifo port=xss

	// Evaluates fitness and local detection
    loop_pso_evalfit: for(unsigned int i = 0; i < n_S; i++) {
        // std::cout << std::endl;
// #pragma UNROLL
		System<_pso_hw_real,_Nh, _Nx, _n_U, _Nu> 
		current_system(
			_u_max, 
			_u_min, 
			_du_max,
			_controlled_state,
			_state_upper_limits, 
			_state_lower_limits, 
			_Q, 
			_Qf, 
			_R, 
			_uss, 
			_Ts
			// ,
			// my_model
		);
// #pragma HLS unroll factor=5 skip_exit_check
        //if(valid_particle[i] == 1){
		// current_system->nmpc_cost_function_topflow(x_curr_local, &x[i][0], xref, &fx[i]);
		current_system.nmpc_cost_function_topflow(local_xref, &local_x[i*part_S], local_xref, &local_fx[i]);
		// fx[i] = rand_real();
		if (local_fx[i] < local_f_ind[i]) {
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nu*_pso_n_U>(&local_y[i*part_S], &local_x[i*part_S]);
			local_f_ind[i] = local_fx[i];           
		}
        //}
    }
    // std::cout << std::endl;

}

// ---------------------------------------------------
#ifdef PSO_CANON

void pso_solver::updateParticles(
//	_pso_hw_real global_min[],

    // _pso_hw_real _x[][_pso_Nu*_pso_n_U],
    // _pso_hw_real _y[][_pso_Nu*_pso_n_U],
    // _pso_hw_real _v[][_pso_Nu*_pso_n_U]
	// _pso_hw_real **_x,
	// _pso_hw_real **_y,
	// _pso_hw_real **_v
){    

    // Update Particles
    for (unsigned int i = 0; i < n_S; ++i) {
		for (unsigned int k = 0; k < Nu; ++k) {
	 	   	for (unsigned int j = 0; j < n_U; ++j) {
    	 		int idx = k*Nu+j;
	
				_pso_hw_real v_tmp = v[i*part_S + idx];
				_pso_hw_real x_tmp = x[i*part_S + idx];
				
				_pso_hw_real r1, r2;
				rand_real(r1); //random->read();
				rand_real(r2); //random->read();

				// v = w*v + c1*r1*(y-x) + c2*r2*(global_min - x)
				_pso_hw_real v_new = w*v_tmp + c1*r1*(y[i*part_S + idx]-x_tmp) + c2*r2*(global_min[idx] - x_tmp);

                v_new = (v_new > max_v) ? max_v : v_new;
				v_new = (v_new < min_v) ? min_v : v_new;
				v[i*part_S + idx] = v_new;

				_pso_hw_real x_new = x_tmp + v_new;
				x[i*part_S + idx] = verifyControlConstrains(x_new, k);
            } // END FOR j
        } // End FOR k
    } // END FOR i
}
#endif

// ---------------------------------------------------

void pso_solver::updateParticlesWithDuConstrains(
	_pso_hw_real *local_u_curr,

	_pso_hw_real local_w,

	_pso_hw_real *local_x,
	_pso_hw_real *local_y,
	_pso_hw_real *local_v,

	_pso_hw_real *local_global_min,

	_pso_hw_real *local_x_max_first,
	_pso_hw_real *local_x_min_first,

	_pso_hw_real *local_du_max,
	_pso_hw_real *local_du_min,
	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min
){

	// Update particles
	_pso_hw_real x_ant[_pso_n_U], x_ant_tmp[_pso_n_U];
    for (unsigned int i = 0; i < n_S; ++i) {
// #pragma UNROLL
// #pragma HLS ALLOCATION operation instances=hmul 	  limit=1 
// #pragma HLS ALLOCATION operation instances=hadd 	  limit=1 
// #pragma HLS ALLOCATION operation instances=hsub 	  limit=1 
// #pragma HLS ALLOCATION function  instances=rand_real limit=1 
		memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant_tmp, (_pso_hw_real*)local_u_curr);
        //if(valid_particle[i] == 1){
    	 for (unsigned int j = 0; j < Nu; ++j) {
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant, (_pso_hw_real*)x_ant_tmp);
    	 	for (unsigned int k = 0; k < n_U; ++k){
// #pragma HLS pipeline II=11 rewind
    	 		// int idx = k*Nu+j;
				int idx = j*n_U+k;
                _pso_hw_real r1,r2;
				pso_solver::rand_real(r1); //random->read();
                pso_solver::rand_real(r2); //random->read();

	            _pso_hw_real v_tmp = local_v[i*part_S + idx];
	            _pso_hw_real x_tmp = local_x[i*part_S + idx];
	            // v = w*v + c1*r1*(y-x) + c2*r2*(global_min - x)
				_pso_hw_real v_new = local_w*v_tmp + _pso_c1*r1*(local_y[i*part_S + idx]-x_tmp) + \
										_pso_c2*r2*(local_global_min[idx] - x_tmp);
	            v_new = (v_new > _pso_max_v) ? _pso_max_v : v_new;
	            v_new = (v_new < _pso_min_v) ? _pso_min_v : v_new;
				
				_pso_hw_real x_ant_k = x_ant[k];
				_pso_hw_real x_new = x_tmp + v_new;
				// _pso_hw_real x_dx_tmp = (j==0) ? x_new : x_new - x_ant_k; //x[i][idx]-x[i][idx-1];
				_pso_hw_real x_dx_tmp = x_new - x_ant_k; //x[i][idx]-x[i][idx-1];
				_pso_hw_real cmp_max = 	(j==0) ? local_x_max_first[k] : local_du_max[k];
				_pso_hw_real cmp_min = 	(j==0) ? local_x_min_first[k] : local_du_min[k];
				x_new = (x_dx_tmp > cmp_max) ? x_ant_k + cmp_max : x_new;
				x_new = (x_dx_tmp < cmp_min) ? x_ant_k + cmp_min : x_new;
                x_new = pso_solver::verifyControlConstrains(
					local_u_max,
					local_u_min,
					x_new, 
					k
				);
				
				local_x[i*part_S + idx] = x_new;
				local_v[i*part_S + idx] = v_new;

				x_ant_tmp[k] = x_new;
	            
            } // END FOR k
        } // END FOR j
       //} // END IF
    } // END FOR i

}
