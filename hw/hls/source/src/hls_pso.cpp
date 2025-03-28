#include "hls_pso.hpp"

#include "hls_system.hpp"
#include "hls_pseudorand.hpp"
#include "aux_functions.hpp"

// using namespace pso_solver;

const unsigned int n_S = _pso_n_S; 	// number of particles
//System * controled_system;

// const int kpso = _KPSO;
// const int stable_zero;
const int particle_last_best = 0;
const int particle_stable_zero = 1;

// PSO Parameters
const unsigned int	pso_maxiter = _pso_maxiter;
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
const unsigned int size_xref = _pso_Nx*_pso_Nh;
const unsigned int part_S = Nu*n_U;
const unsigned int part_S_mem = _pso_n_S * _pso_Nu*_pso_n_U;
// _pso_hw_real u_from_parameters[_N*_pso_n_U];

// Pseudo Random Core Generator
typedef pseudoRand_gen<_pso_hw_real> _pso_randCore_t;
// _pso_randCore_t randGen(
// 	(const float)pso_rand_min, 
// 	(const float)pso_rand_max
// );
_pso_randCore_t rand_core;

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
){

	// _pso_hw_real w;

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
	_pso_hw_real bestfitness[_pso_maxiter];
	_pso_hw_real global_min[_pso_Nu*_pso_n_U];
	// _pso_hw_real _pso_x_min_first[_pso_n_U], _pso_x_max_first[_pso_n_U];
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
	#ifdef __SYNTHESIS__
		_rand_real_stream __rand_port;
		// _pso_hw_real rand_value;
		// rand_real(rand_value);
		// __rand_port.write(rand_value)
		memcpy_loop_rolled<_pso_hw_real, volatile _pso_hw_real, _pso_n_U>(u_curr_local, 	u_curr);
		memcpy_loop_rolled<_pso_hw_real, volatile _pso_hw_real, _pso_Nx>(x_curr_local, 	x_curr);
		memcpy_loop_rolled<_pso_hw_real, volatile _pso_hw_real, _pso_Nx*_pso_Nh>(xref_local, xref);
	#endif

	initializeParticles_set(
#ifndef __SYNTHESIS__
		u_curr,
		x_curr,
		xref,
#endif
		last_best,
		
		u_curr_local,
		x_curr_local,
		xref_local,

		x,y,v,

		f_ind,

		_pso_du_min,
		_pso_du_max,
		(_pso_hw_real *)_pso_u_min,
		(_pso_hw_real *)_pso_u_max,
		(_pso_hw_real *)_pso_uss_local

#ifdef __SYNTHESIS__
		, __rand_port
#endif
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
        evaluateFitnessAndDetectLocalBest(
			x, 
			y, 
			
			x_curr_local, 

			xref_local,
			f_ind
			// ,fx
		);

        // Global Minimum detection
		// best_pos = detectGlobalMinimum(k);
		detectGlobalMinimum(
			f_ind,
			bestfitness,
			global_min,
			k,
			y
		);
		// memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nu*_pso_n_U>(global_min, (_pso_hw_real *)&y[best_pos * part_S]);

#ifndef PSO_CANON
        updateParticlesWithDuConstrains(
			u_curr_local,
			// w,

			x,
			y,
			v,

			global_min,

			// _pso_x_max_first,
			// _pso_x_min_first,

			_pso_du_max,
			_pso_du_min,
			(_pso_hw_real *)_pso_u_max,
			(_pso_hw_real *)_pso_u_min
#ifdef __SYNTHESIS__
			, __rand_port
#endif
		);    
#else
		updateParticles(
			x,
			y,
			v,

			global_min,

			(_pso_hw_real *)_pso_u_max,
			(_pso_hw_real *)_pso_u_min
		);
#endif
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
		
		// w+= _pso_slope;
	}


	// Return best value
    memcpy_loop_rolled<_pso_hw_real,_pso_hw_real,_pso_Nu*_pso_n_U>(new_best, global_min);
	// for (unsigned int i = 0; i < Nu*n_U; ++i) {
    //     new_best[i] = global_min[i];
    // }
    
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

// ---------------------------------------------------

int pso_fsm(
	// Top Level 
	_pso_hw_real u_curr[_pso_n_U],
	_pso_hw_real x_curr[_pso_Nx],
	_pso_hw_real xref[_pso_Nx*_pso_Nh],
	_pso_hw_real last_best[_pso_Nu*_pso_n_U],

	_pso_hw_real new_best[_pso_Nu*_pso_n_U],
	_pso_hw_real bestfitness_top[_pso_maxiter],

	// Local memories for system constraints created now
	_pso_hw_real last_best_local[_pso_Nu*_pso_n_U],
	_pso_hw_real u_curr_local[_pso_n_U],
	_pso_hw_real x_curr_local[_pso_Nx],
	_pso_hw_real xref_local[_pso_Nx*_pso_Nh],

	// State Machine signals
	bool &rand_core_ap_start,
	bool &init_s_ap_start,
	bool init_s_ap_idle,
	bool &eval_s_ap_start,
	bool eval_s_ap_idle,
	bool &detect_min_ap_start,
	bool detect_min_ap_idle,
	bool &update_s_ap_start,
	bool update_s_ap_idle,

	bool &rst_cores,

	// Particle Variables
	uint32_t &k,
	_pso_hw_real local_bestfitness[_pso_maxiter],
	_pso_hw_real global_min[_pso_Nu*_pso_n_U]
	
){
#pragma HLS INTERFACE mode=s_axilite port=return bundle=control

#pragma HLS INTERFACE mode=ap_none port=rand_core_ap_start
#pragma HLS INTERFACE mode=ap_none port=init_s_ap_start
#pragma HLS INTERFACE mode=ap_none port=init_s_ap_idle
#pragma HLS INTERFACE mode=ap_none port=eval_s_ap_start
#pragma HLS INTERFACE mode=ap_none port=eval_s_ap_idle
#pragma HLS INTERFACE mode=ap_none port=detect_min_ap_start
#pragma HLS INTERFACE mode=ap_none port=detect_min_ap_idle
#pragma HLS INTERFACE mode=ap_none port=update_s_ap_start
#pragma HLS INTERFACE mode=ap_none port=update_s_ap_idle
#pragma HLS INTERFACE mode=ap_none port=rst_cores

#pragma HLS interface mode=ap_vld register port=k

#pragma HLS INTERFACE mode=m_axi port=u_curr			offset=slave 	depth=n_U			bundle=nmpc_io
#pragma HLS INTERFACE mode=m_axi port=x_curr			offset=slave 	depth=Nx			bundle=nmpc_io 
#pragma HLS INTERFACE mode=m_axi port=xref				offset=slave 	depth=size_xref		bundle=nmpc_io
#pragma HLS INTERFACE mode=m_axi port=last_best			offset=slave 	depth=part_S		bundle=nmpc_io

#pragma HLS INTERFACE mode=m_axi port=new_best			offset=slave 	depth=part_S		bundle=nmpc_io
#pragma HLS INTERFACE mode=m_axi port=bestfitness_top	offset=slave 	depth=pso_maxiter	bundle=nmpc_io

#pragma HLS INTERFACE mode=m_axi port=u_curr_local		offset=direct	depth=n_U			bundle=current_mem 
#pragma HLS INTERFACE mode=m_axi port=x_curr_local		offset=direct	depth=Nx			bundle=current_mem 
#pragma HLS INTERFACE mode=m_axi port=xref_local		offset=direct	depth=size_xref		bundle=current_mem 
#pragma HLS INTERFACE mode=m_axi port=last_best_local	offset=direct	depth=part_S		bundle=current_mem 

#pragma HLS INTERFACE mode=m_axi port=local_bestfitness	offset=direct 	depth=pso_maxiter	bundle=pso_mem
#pragma HLS INTERFACE mode=m_axi port=global_min		offset=direct 	depth=part_S		bundle=pso_mem

	uint32_t k_local = 0;
	k = 0;
	uint16_t state = 0;
	
	rst_cores			= false;
	rand_core_ap_start	= false;
	init_s_ap_start		= false;
	eval_s_ap_start 	= false;
	detect_min_ap_start = false;
	update_s_ap_start 	= false;

	bool flag=true;
	do{
		switch (state)
		{
		case 0: // Save Input data in memory banks
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(u_curr_local, 	u_curr);
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nx>(x_curr_local, 	x_curr);
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nx*_pso_Nh>(xref_local, xref);
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nu*_pso_n_U>(last_best_local, last_best);

			rst_cores			= true;
			rand_core_ap_start	= false;
			init_s_ap_start		= false;
			eval_s_ap_start 	= false;
			detect_min_ap_start = false;
			update_s_ap_start 	= false;
						
			state = 1;
			break;

		case 1:
			rst_cores			= false;
			rand_core_ap_start	= true;
			init_s_ap_start		= true;
			eval_s_ap_start 	= false;
			detect_min_ap_start = false;
			update_s_ap_start 	= false;

			state = 2;
			break;

		case 2: // Init State
			rst_cores			= false;
			rand_core_ap_start	= true;
			init_s_ap_start		= false;
			eval_s_ap_start 	= false;
			detect_min_ap_start = false;
			update_s_ap_start 	= false;
			
			if (init_s_ap_idle)
				state = 2;
			else
				state = 3;

			break;

		case 3: // Interactive Start
			rst_cores			= false;
			rand_core_ap_start	= false;
			init_s_ap_start		= false;
			eval_s_ap_start 	= true;
			detect_min_ap_start = false;
			update_s_ap_start 	= false;
			
			state = 4;
			break;

		case 4: // Cost F
			rst_cores			= false;
			rand_core_ap_start	= false;
			init_s_ap_start		= false;
			eval_s_ap_start 	= false;
			detect_min_ap_start = false;
			update_s_ap_start 	= false;
			
			if (eval_s_ap_idle)
				state = 4;
			else
				state = 5;
			break;

		case 5:// Detect Min
			rst_cores			= false;
			rand_core_ap_start	= false;
			init_s_ap_start		= false;
			eval_s_ap_start 	= false;
			detect_min_ap_start = true;
			update_s_ap_start 	= false;
			
			state = 6;
			break;

		case 6:// 
			rst_cores			= false;
			rand_core_ap_start	= false;
			init_s_ap_start		= false;
			eval_s_ap_start 	= false;
			detect_min_ap_start = false;
			update_s_ap_start 	= false;
			
			if (detect_min_ap_idle)
				state = 6;
			else
				state = 7;
			break;

		case 7://  Update S
			rst_cores			= false;
			rand_core_ap_start	= true;
			init_s_ap_start		= false;
			eval_s_ap_start 	= false;
			detect_min_ap_start = false;
			update_s_ap_start 	= true;
			
			state = 8;
			break;
		
		case 8:// 
			rst_cores			= false;
			rand_core_ap_start	= true;
			init_s_ap_start		= false;
			eval_s_ap_start 	= false;
			detect_min_ap_start = false;
			update_s_ap_start 	= false;
			
			if (update_s_ap_idle)
				state = 8;
			else
				k_local++;
				k = k_local;
				if (k_local == _pso_maxiter)
					state = 9;
				else
					state = 3;
			break;

		case 9:// Finish and save data
			rst_cores			= false;
			rand_core_ap_start	= false;
			init_s_ap_start		= false;
			eval_s_ap_start 	= false;
			detect_min_ap_start = false;
			update_s_ap_start 	= false;
			
			state = 10;
			memcpy_loop_rolled<_pso_hw_real,_pso_hw_real,_pso_Nu*_pso_n_U>(new_best, global_min);
			memcpy_loop_rolled<_pso_hw_real,_pso_hw_real,_pso_maxiter>(bestfitness_top,local_bestfitness);
			flag = false;

			break;
		
		case 10:// Finish and save data
			flag = false;
			state = 10;
			break;

		default:
			state = 1;
			break;
		}
	}while(flag);

	return k_local;
}

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
){

#pragma HLS INTERFACE mode=m_axi port=local_x offset=off depth=part_S_mem bundle=x_mem
#pragma HLS INTERFACE mode=m_axi port=local_y offset=off depth=part_S_mem bundle=y_mem
#pragma HLS INTERFACE mode=m_axi port=local_v offset=off depth=part_S_mem bundle=v_mem

#pragma HLS INTERFACE mode=m_axi port=u_curr_local	offset=direct	depth=n_U		bundle=current_mem 
#pragma HLS INTERFACE mode=m_axi port=x_curr_local	offset=direct	depth=Nx		bundle=current_mem 
#pragma HLS INTERFACE mode=m_axi port=xref_local	offset=direct	depth=size_xref	bundle=current_mem 
#pragma HLS INTERFACE mode=m_axi port=last_best 	offset=direct	depth=part_S 	bundle=current_mem 

#pragma HLS INTERFACE mode=m_axi port=f_ind_local 	offset=direct	depth=n_S bundle=pso_mem

#pragma HLS INTERFACE mode=m_axi port=local_du_min 	offset=direct	depth=n_U bundle=constraints_mem
#pragma HLS INTERFACE mode=m_axi port=local_du_max 	offset=direct	depth=n_U bundle=constraints_mem
#pragma HLS INTERFACE mode=m_axi port=local_u_min 	offset=direct	depth=n_U bundle=constraints_mem
#pragma HLS INTERFACE mode=m_axi port=local_u_max 	offset=direct	depth=n_U bundle=constraints_mem
#pragma HLS INTERFACE mode=m_axi port=local_uss 	offset=direct	depth=n_U bundle=constraints_mem

#ifdef __SYNTHESIS__
#pragma HLS interface mode=axis port=__rand_port 
#endif

#ifndef __SYNTHESIS__
	memcpy_loop_rolled<_pso_hw_real, volatile _pso_hw_real, _pso_n_U>(u_curr_local, 	u_curr);
	memcpy_loop_rolled<_pso_hw_real, volatile _pso_hw_real, _pso_Nx>(x_curr_local, 	x_curr);
	memcpy_loop_rolled<_pso_hw_real, volatile _pso_hw_real, _pso_Nx*_pso_Nh>(xref_local, xref);
#endif

#ifndef PSO_CANON
	// w = _pso_w0;
	calculate_du_min(
		local_du_min,
		local_du_max
	);

    initializeParticlesWithDuConstrains(
		u_curr_local,
		local_du_max,
		
		local_u_max,
		local_u_min,

		local_x,
		local_y,
		local_v
		
		// ,valid_particle

#ifdef __SYNTHESIS__
		,__rand_port
#endif		
	);

	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real,(Nu*n_U)>(&local_x[particle_last_best*part_S], last_best);

	initializeStableZero(
		local_uss,
		u_curr_local, 

		local_u_max,
		local_u_min,
		local_du_max,
		local_du_min,
		
		local_x,

		particle_stable_zero
	);
#else
	initializeParticles(
		local_u_max,
		local_u_min,

		local_x,
		local_y,
		local_v
	);
		
#endif

    // initializeBestLocalFitness();
	memset_loop<_pso_hw_real>(f_ind_local, (const _pso_hw_real)H_MAX, n_S);
}
// ---------------------------------------------------

void rand_real(_pso_hw_real &rand_value){
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

void  detectGlobalMinimum(
	_pso_hw_real local_find[_pso_n_S], 
	_pso_hw_real local_bestfitness[_pso_maxiter],
	_pso_hw_real local_global_min[_pso_Nu*_pso_n_U],
	int interaction,
	
	_pso_hw_real local_y[_pso_n_S * _pso_Nu*_pso_n_U]
){
#pragma HLS INTERFACE mode=m_axi port=local_bestfitness	offset=direct bundle=current_mem depth=pso_maxiter
#pragma HLS INTERFACE mode=m_axi port=local_global_min	offset=direct bundle=current_mem depth=n_U

#pragma HLS INTERFACE mode=m_axi port=local_find 	offset=direct	depth=n_S bundle=constraints_mem

#pragma HLS INTERFACE mode=m_axi port=local_y offset=off depth=part_S_mem bundle=y_mem

#pragma HLS interface mode=ap_vld register port=interaction

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

// void initializeConstrains(
// 	volatile _pso_hw_real *u_curr,
// 	_pso_hw_real *local_du_max,
// 	_pso_hw_real *local_du_min,
// 	_pso_hw_real *local_u_max,
// 	_pso_hw_real *local_u_min,

// 	_pso_hw_real *local_x_max_first,
// 	_pso_hw_real *local_x_min_first
// ){
// 	// Initialize constrains based on current statecmp
// #pragma HLS inline
// 	initializeConstrains_loop: for (unsigned int i = 0; i < n_U; ++i){
// #pragma HLS pipeline
// 		_pso_hw_real u_curr_tmp = u_curr[i];
// 		_pso_hw_real local_du_max_i = local_du_max[i];
// 		_pso_hw_real local_du_min_i = local_du_min[i];

//         _pso_hw_real x_max_first_temp = u_curr_tmp + local_du_max_i; //controled_system->acc_max[i];
// 		_pso_hw_real x_min_first_temp = u_curr_tmp + local_du_min[i];//controled_system->acc_min[i];        

// 		if (x_max_first_temp > local_u_max[i])
// 			local_x_max_first[i] = ;
// 		else
// 			local_x_max_first[i] = local_u_max[i];

// 		if (x_min_first_temp < local_u_min[i])
// 			local_x_min_first[i] = ;
// 		else
// 			local_x_min_first[i] = local_u_min[i];	
// 	}
// }

// ---------------------------------------------------

void calculate_du_min(
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

void equalizeParticles(
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

_pso_hw_real verifyControlConstrains(
	_pso_hw_real local_u_max,
	_pso_hw_real local_u_min,
	_pso_hw_real value
){
#pragma HLS inline
	_pso_hw_real return_value = value;

	if (return_value > local_u_max)
		return_value = local_u_max;
	else if (return_value < local_u_min)
		return_value = local_u_min;
	else
		return_value = return_value;

	return return_value;
}

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
){
// INITIALIZATION OF PARTICLES WITH Delta u CONTRAINS FOR MPC
#pragma HLS inline
// #pragma HLS ALLOCATION instances=rand_real limit=1 function
// #pragma HLS ALLOCATION instances=hmul limit=1 operation
// #pragma HLS ALLOCATION instances=hadd limit=1 operation
// #pragma HLS ALLOCATION instances=hsub limit=1 operation

	_pso_hw_real x_ant[_pso_n_U];
	_pso_hw_real x_ant_tmp[_pso_n_U];

	initializeS_du_loop_S:for (unsigned int i = 0; i < n_S; ++i) {
		memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant_tmp, (_pso_hw_real *)u_curr);
		initializeS_du_loop_N:for (unsigned int j = 0; j < Nu; ++j) {
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant, (_pso_hw_real *)x_ant_tmp);
			initializeS_du_loop_x:for (unsigned int k = 0; k < n_U; ++k) {
#pragma HLS pipeline
				int idx = j*n_U + k;
				_pso_hw_real rand_tmp;
				_pso_hw_real local_du_max_k = local_du_max[k];
#ifdef __SYNTHESIS__
				rand_tmp = __rand_port.read();
#else
				rand_real(rand_tmp);
#endif			
				_pso_hw_real du_tmp = rand_tmp*(_pso_hw_real)2.0*local_du_max_k - local_du_max_k;
		        _pso_hw_real x_tmp = x_ant[k] + du_tmp; //random->read(); 
		        // _pso_hw_real x_tmp = ( x_ant[k] + (-du_max[k]) ) + ((_pso_hw_real)2.0*du_max[k]) * rand_tmp; //random->read(); 
				x_tmp = verifyControlConstrains(
					local_u_max[k],
					local_u_min[k],
					x_tmp
				);
                local_x[i*part_S + idx] = x_tmp;
				local_y[i*part_S + idx] = x_tmp; // uss_local[k];
		        local_v[i*part_S + idx] = _pso_init_v;
                // valid_particle[i] = true;
				x_ant_tmp[k] = x_tmp;
				idx++;
		    }
		}
	}
}

// ---------------------------------------------------

void  initializeParticles(
	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min,

	_pso_hw_real *local_x,
	_pso_hw_real *local_y,
	_pso_hw_real *local_v
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
				_pso_hw_real x_tmp;
				if (j == 0) 
					x_tmp = local_u_min[k] + (local_u_max[k]-local_u_min[k]) * rand_tmp; //random->read();
				else
					x_tmp = x_ant[k] - d_max + (d2_max) * rand_tmp; //random->read();
				local_x[i*part_S + idx] = x_tmp; 
				local_y[i*part_S + idx] = x_tmp;
				local_v[i*part_S + idx] = _pso_init_v;  
				x_ant[k] = x_tmp;
			}
        }
        // valid_particle[i] = 1;
    }
}


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
#pragma HLS pipeline
        	idx = k*n_U + i;
		    _pso_hw_real x_tmp = uss_local[i];
			// _pso_hw_real comp_tmp = (k == 0) ? x_tmp - u_curr[i] : x_tmp - x_ant[i];
			_pso_hw_real comp_tmp = x_tmp - x_ant[i] ;
			
			if (comp_tmp > local_du_max[i])
				x_tmp = local_du_max[i];
			else if (comp_tmp < local_du_min[i])
				x_tmp = local_du_min[i];
			else
				x_tmp = x_tmp;

            local_x[index*part_S + idx] = verifyControlConstrains(
				local_u_max[k],
				local_u_min[k],
				x_tmp
			);
			x_ant_tmp[i] = x_tmp;
			idx++;
        }
   	}
}

// ---------------------------------------------------

void evaluateFitnessAndDetectLocalBest(
	_pso_hw_real local_x[_pso_n_S * _pso_Nu*_pso_n_U],
	_pso_hw_real local_y[_pso_n_S * _pso_Nu*_pso_n_U],

	_pso_hw_real local_x_curr[_pso_Nx],

	_pso_hw_real local_xref[_pso_Nx*_pso_Nu], 
	_pso_hw_real local_f_ind[_pso_n_U]
	// _pso_hw_real *local_fx//[_pso_n_U],
){

#pragma HLS INTERFACE mode=m_axi port=local_x offset=off depth=part_S_mem bundle=x_mem
#pragma HLS INTERFACE mode=m_axi port=local_y offset=off depth=part_S_mem bundle=y_mem

#pragma HLS INTERFACE mode=m_axi port=local_x_curr	offset=direct	depth=Nx		bundle=current_mem 
#pragma HLS INTERFACE mode=m_axi port=local_xref	offset=direct	depth=size_xref	bundle=current_mem 
#pragma HLS INTERFACE mode=m_axi port=local_f_ind 	offset=direct	depth=n_S bundle=constraints_mem

	_pso_hw_real fx[_pso_n_S];
#pragma HLS array_partition variable=fx	type=complete factor=n_S 

	// Evaluates fitness and local detection
    loop_pso_evalfit: for(unsigned int i = 0; i < n_S; i++) {
		System<_pso_hw_real, _hw_model_real, _Nh, _Nx, _n_U, _Nu> 
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
		current_system.nmpc_cost_function_topflow(local_x_curr, &local_x[i*part_S], local_xref, &fx[i]);
	}
	
	for (int i = 0; i < n_S; i++){	
#pragma UNROLL off
		if (fx[i] < local_f_ind[i]) {
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nu*_pso_n_U>(&local_y[i*part_S], &local_x[i*part_S]);
			local_f_ind[i] = fx[i] ;
		}
    }
}

// ---------------------------------------------------

void updateParticles(
	_pso_hw_real *local_x,
	_pso_hw_real *local_y,
	_pso_hw_real *local_v,

	_pso_hw_real *local_global_min,

	_pso_hw_real *local_u_max,
	_pso_hw_real *local_u_min
){    
	static _pso_hw_real w = _pso_w0;

    // Update Particles
    for (unsigned int i = 0; i < n_S; ++i) {
		for (unsigned int j = 0; j < Nu; ++j) {
	 	   	for (unsigned int k = 0; k < n_U; ++k) {
    	 		int idx = j*n_U+k;
	
				_pso_hw_real v_tmp = local_v[i*part_S + idx];
				_pso_hw_real x_tmp = local_x[i*part_S + idx];
				
				_pso_hw_real r1, r2;
				rand_real(r1); //random->read();
				rand_real(r2); //random->read();

				// v = w*v + c1*r1*(y-x) + c2*r2*(global_min - x)
				_pso_hw_real v_new = w*v_tmp + _pso_c1*r1*(local_y[i*part_S + idx]-x_tmp) +\
										_pso_c2*r2*(local_global_min[idx] - x_tmp);

                if (v_new > _pso_max_v) 
					v_new = _pso_max_v;
				else if (v_new < _pso_min_v)
					v_new = _pso_min_v;
				else 
					v_new = v_new;
				local_v[i*part_S + idx] = v_new;

				_pso_hw_real x_new = x_tmp + v_new;
				local_x[i*part_S + idx] = verifyControlConstrains(
					local_u_max[k],
					local_u_min[k],
					x_new);
            } // END FOR j
        } // End FOR k
    } // END FOR i

	_pso_hw_real w_tmp = w + _pso_slope;
	if (w_tmp <= _pso_wf)
		w = _pso_w0;
	else
		w = w_tmp;

	return;
}

// ---------------------------------------------------

void updateParticlesWithDuConstrains(
	_pso_hw_real local_u_curr[_pso_n_U],

	// _pso_hw_real local_w,

	_pso_hw_real local_x[_pso_n_S * _pso_Nu*_pso_n_U],
	_pso_hw_real local_y[_pso_n_S * _pso_Nu*_pso_n_U],
	_pso_hw_real local_v[_pso_n_S * _pso_Nu*_pso_n_U],

	_pso_hw_real local_global_min[_pso_Nu*_pso_n_U],

	// _pso_hw_real *local_x_max_first,
	// _pso_hw_real *local_x_min_first,

	_pso_hw_real local_du_max[_pso_n_U],
	_pso_hw_real local_du_min[_pso_n_U],
	_pso_hw_real local_u_max[_pso_n_U],
	_pso_hw_real local_u_min[_pso_n_U]

#ifdef __SYNTHESIS__
	,_rand_real_stream &__rand_port
#endif
){
#pragma HLS INTERFACE mode=m_axi port=local_x offset=off depth=part_S_mem bundle=x_mem
#pragma HLS INTERFACE mode=m_axi port=local_y offset=off depth=part_S_mem bundle=y_mem
#pragma HLS INTERFACE mode=m_axi port=local_v offset=off depth=part_S_mem bundle=v_mem

#pragma HLS INTERFACE mode=m_axi port=local_u_curr		offset=direct	depth=n_U		bundle=current_mem 
#pragma HLS INTERFACE mode=m_axi port=local_global_min	offset=direct	depth=n_U		bundle=current_mem

#pragma HLS INTERFACE mode=m_axi port=local_du_min 		offset=direct	depth=n_U bundle=constraints_mem
#pragma HLS INTERFACE mode=m_axi port=local_du_max 		offset=direct	depth=n_U bundle=constraints_mem
#pragma HLS INTERFACE mode=m_axi port=local_u_min 		offset=direct	depth=n_U bundle=constraints_mem
#pragma HLS INTERFACE mode=m_axi port=local_u_max 		offset=direct	depth=n_U bundle=constraints_mem

#ifdef __SYNTHESIS__
#pragma HLS interface mode=axis port=__rand_port 
#endif
	static _pso_hw_real w = _pso_w0;
	
	_pso_hw_real __local_global_min[_pso_Nu*_pso_n_U];
	_pso_hw_real __local_du_max[_pso_n_U];
	_pso_hw_real __local_du_min[_pso_n_U];
	_pso_hw_real __local_u_max[_pso_n_U];
	_pso_hw_real __local_u_min[_pso_n_U];

	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nu*_pso_n_U>(__local_global_min, local_global_min);
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(__local_du_max, local_du_max);
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(__local_du_min, local_du_min);
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(__local_u_max, local_u_max);
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(__local_u_min, local_u_min);

	// Update particles
	_pso_hw_real x_ant[_pso_n_U], x_ant_tmp[_pso_n_U];
    for (unsigned int i = 0; i < n_S; ++i) {	
		memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant_tmp, (_pso_hw_real*)local_u_curr);
        //if(valid_particle[i] == 1){
    	 for (unsigned int j = 0; j < Nu; ++j) {
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant, (_pso_hw_real*)x_ant_tmp);
    	 	for (unsigned int k = 0; k < n_U; ++k){
#pragma HLS pipeline II=3 rewind
    	 		// int idx = k*Nu+j;
				int idx = j*n_U+k;
                _pso_hw_real r1,r2;
#ifdef __SYNTHESIS__
				r1 = __rand_port.read();
				r2 = __rand_port.read();
#else
				rand_real(r1); //random->read();
                rand_real(r2); //random->read();
#endif		
	            _pso_hw_real v_tmp = local_v[i*part_S + idx];
	            _pso_hw_real x_tmp = local_x[i*part_S + idx];
	            // v = w*v + c1*r1*(y-x) + c2*r2*(global_min - x)
				_pso_hw_real v_new = w*v_tmp + _pso_c1*r1*(local_y[i*part_S + idx] - x_tmp) + \
										_pso_c2*r2*(__local_global_min[idx] - x_tmp);
	            if (v_new > _pso_max_v) 
					v_new = _pso_max_v;
				else if (v_new < _pso_min_v) 
					v_new = _pso_min_v;
				else 
					v_new = v_new;
				
				_pso_hw_real x_ant_k = x_ant[k];
				_pso_hw_real x_new = x_tmp + v_new;
				// _pso_hw_real x_dx_tmp = (j==0) ? x_new : x_new - x_ant_k; //x[i][idx]-x[i][idx-1];
				_pso_hw_real x_dx_tmp = x_new - x_ant_k; //x[i][idx]-x[i][idx-1];
				_pso_hw_real cmp_max = __local_du_max[k];
				_pso_hw_real cmp_min = __local_du_min[k];

				if (x_dx_tmp > cmp_max)
					x_new = x_ant_k + cmp_max;
				else if (x_dx_tmp < cmp_min)
					x_new = x_ant_k + cmp_min;
				else
					x_new = x_new;

                x_new = verifyControlConstrains(
					__local_u_max[k],
					__local_u_min[k],
					x_new
				);
				
				local_x[i*part_S + idx] = x_new;
				local_v[i*part_S + idx] = v_new;

				x_ant_tmp[k] = x_new;
	            
            } // END FOR k
        } // END FOR j
       //} // END IF
    } // END FOR i

	_pso_hw_real w_tmp = w + _pso_slope;
	if (w_tmp <= _pso_wf)
		w = _pso_w0;
	else
		w = w_tmp;

	return;
}
