#ifndef HLS_PSO_HPP
#define HLS_PSO_HPP

//#define DEBUG_HLS

//#include "config.hpp"
#include "hls_system.hpp"
#include "hls_pseudorand.hpp"
#include "aux_functions.hpp"

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

template<
    class _pso_hw_real,
    // class _pso_randCore_t,
    // class _pso_system_t,
	unsigned _pso_n_S,
	unsigned _pso_maxiter,
    unsigned _pso_Nh,
    unsigned _pso_Nx,
    unsigned _pso_n_U,
    unsigned _pso_Nu
>class PSO{
protected:
	// Pseudo random generator variables
	

	const unsigned int n_S = _pso_n_S; 	// number of particles
	//System * controled_system;

	// const int kpso = _KPSO;
	const int stable_zero;
	const uint8_t particle_last_best = 0;
	const uint8_t particle_stable_zero = 1;

	// PSO Parameters
	const unsigned int	maxiter = _pso_maxiter;
	const _pso_hw_real 	max_v;
	const _pso_hw_real 	min_v;
	const _pso_hw_real 	w0; // initial weight
	const _pso_hw_real	wf; // final weight
	const _pso_hw_real	slope;
	const _pso_hw_real	c1; // cognitive coefficient
	const _pso_hw_real	c2; // social coefficient
	
	_pso_hw_real w;
	
	_pso_hw_real threshold;
	int stop_criteria;

	const _pso_hw_real *u_min;//[_pso_n_U];
	const _pso_hw_real *u_max;//[_pso_n_U];
	const _pso_hw_real *du_max;//[_pso_n_U];
	_pso_hw_real du_min[_pso_n_U];//[_pso_n_U];
	const _pso_hw_real *uss_local;//[_pso_n_U];
	_pso_hw_real x_max_first[_pso_n_U];
	_pso_hw_real x_min_first[_pso_n_U];

	const _pso_hw_real init_v; // initial velocity

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
	int valid_particle[_pso_n_S];
	int number_of_active_particles;

	_pso_hw_real f_ind[_pso_n_S];
	_pso_hw_real fx[_pso_n_S];
	_pso_hw_real bestfitness[_pso_maxiter];
	_pso_hw_real global_min[_pso_Nu*_pso_n_U];
	//_pso_hw_real previous_control[_Parametrization];

	// Controlled System Configuration
	const unsigned int Nh = _pso_Nh;
	const unsigned int Nx = _pso_Nx;
	const unsigned int n_U = _pso_n_U;
	const unsigned int Nu = _pso_Nu;
	const unsigned int part_S = Nu*n_U;
	// _pso_hw_real u_from_parameters[_N*_pso_n_U];

	// _pso_system_t current_system;
	// _pso_randCore_t randGen;

	const float pso_rand_min = 0.0f;
	const float pso_rand_max = 1.0f;
	
	// Pseudo Random Core Generator
	typedef pseudoRand_gen _pso_randCore_t;
public:
// Init Execute
constexpr PSO(
	// PSO Configuration
	const int 	 	__stable_zero,
	const _pso_hw_real  __max_v,
	const _pso_hw_real  __w0,
	const _pso_hw_real 	__wf,
	const _pso_hw_real 	__slope,
	const _pso_hw_real 	__c1,
	const _pso_hw_real 	__c2,
	const _pso_hw_real	*__u_min,
	const _pso_hw_real	*__u_max,
	const _pso_hw_real	*__du_max,
	const _pso_hw_real 	*__uss
/*	
	,
	const unsigned short _controlled_state[_pso_Nx],
	const _pso_hw_real _state_upper_limits[_pso_Nx],
	const _pso_hw_real _state_lower_limits[_pso_Nx],
	const _pso_hw_real _Q[_pso_Nx],
	const _pso_hw_real _Qf[_pso_Nx],
	const _pso_hw_real _R[_pso_n_U],
	const _pso_hw_real _Ts,
	const int _rand_seed[_pso_n_S]
*/
	// ,
	// _pso_system_t __current_system,
	// _pso_randCore_t __randGen
	) : max_v(__max_v), min_v(-__max_v), w0(__w0), wf(__wf), 
		slope(__slope), c1(__c1), c2(__c2),
		stable_zero(__stable_zero), init_v(__max_v*0.1),
		u_min(__u_min), u_max(__u_max), du_max(__du_max), uss_local(__uss)
		// ,
/*		
		Ts(_Ts), controlled_state(_controlled_state), 
        state_upper_limits(_state_upper_limits), state_lower_limits(_state_lower_limits),
        Q(_Q), Qf(_Qf), R(_R), 
*/
		// randGen(__randGen)
		// , current_system(__current_system)
{
}

// ---------------------------------------------------
int execute(
	volatile _pso_hw_real x_curr[_pso_Nx], 
	volatile _pso_hw_real u_curr[_pso_n_U], 
	int iteration, 
	volatile _pso_hw_real last_best[_pso_n_U*_pso_Nu], 
	volatile _pso_hw_real xref[_pso_Nx*_pso_Nh],
	volatile _pso_hw_real uref[_pso_n_U], 
	volatile _pso_hw_real xss[_pso_Nx],
	//volatile _pso_hw_real uss[_pso_n_U], 
	
	_pso_hw_real new_best[_pso_Nu*_pso_n_U],
	_pso_hw_real * J
){

#pragma HLS ALLOCATION operation instances=hadd limit=2
#pragma HLS ALLOCATION operation instances=hsub limit=2
#pragma HLS ALLOCATION operation instances=hmul limit=2

// #pragma HLS bind_storage variable=x_max_first type=RAM_2P impl=LUTRAM
// #pragma HLS bind_storage variable=x_min_first type=RAM_2P impl=LUTRAM

// #pragma HLS bind_storage variable=x type=RAM_T2P impl=BRAM
// #pragma HLS stream variable=x type=shared
// #pragma HLS bind_storage variable=y type=RAM_T2P impl=BRAM
// #pragma HLS stream variable=y type=shared
// #pragma HLS bind_storage variable=v type=RAM_T2P impl=BRAM
// #pragma HLS stream variable=v type=shared

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
	_pso_hw_real u_curr_local[n_U];
#pragma HLS bind_storage variable=u_curr_local type=RAM_2P impl=LUTRAM
// #pragma HLS stream variable=u_curr_local type=shared

	_pso_hw_real x_curr_local[Nx];
#pragma HLS bind_storage variable=x_curr_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=x_curr_local type=shared

	_pso_hw_real xref_local[Nx*Nh];
#pragma HLS bind_storage variable=xref_local type=RAM_2P impl=BRAM
#pragma HLS stream variable=xref_local type=shared

	_pso_hw_real uref_local[n_U];
#pragma HLS bind_storage variable=uref_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=uref_local type=shared

	_pso_hw_real xss_local[Nx];
#pragma HLS bind_storage variable=xss_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=xss_local type=shared
#else
	_pso_hw_real *u_curr_local;
	_pso_hw_real *x_curr_local;
	_pso_hw_real *xref_local;
	_pso_hw_real *uref_local;
	_pso_hw_real *xss_local;
	x_curr_local = (_pso_hw_real *) malloc(Nx*sizeof(_pso_hw_real));
	u_curr_local = (_pso_hw_real *) malloc(n_U*sizeof(_pso_hw_real));
	xref_local = (_pso_hw_real *) malloc(Nx*Nh*sizeof(_pso_hw_real));
	uref_local = (_pso_hw_real *) malloc(n_U*sizeof(_pso_hw_real));
	xss_local = (_pso_hw_real *) malloc(Nx*sizeof(_pso_hw_real));
#endif
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(u_curr_local, 	(_pso_hw_real *)u_curr);
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nx>(x_curr_local, 	(_pso_hw_real *)x_curr);
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nx*_pso_Nh>(xref_local, (_pso_hw_real *)xref);
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(uref_local, 	(_pso_hw_real *)uref);
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nx>(xss_local, 		(_pso_hw_real *)xss);

	w = w0;
	calculate_du_min();

	initializeConstrains(u_curr_local);
    initializeParticlesWithDuConstrains(u_curr_local);

	// equalizeParticles(iteration);

    initializeLastBestKPSO(last_best, particle_last_best);

	initializeStableZero(u_curr_local, particle_stable_zero);

    initializeBestLocalFitness();

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

	while (k < maxiter) {
        evaluateFitnessAndDetectLocalBest(
			// x, 
			// y, 
			x_curr_local, 
			xref_local, 
			uref_local, 
			xss_local
		);

        // Global Minimum detection
		best_pos = detectGlobalMinimum(k);
		memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nu*_pso_n_U>(global_min, (_pso_hw_real *)&y[best_pos * part_S]);

        updateParticlesWithDuConstrains();    
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
		w+= slope;
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
	if(uref_local != NULL) free(uref_local);
	if(xss_local != NULL) free(xss_local);
#endif
	return k;
}

private:
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


_pso_hw_real rand_real(){
#if (defined(__SYNTHESIS__) || defined(SYNTH_RAND))
	_pso_randCore_t randGen(
		(const float)pso_rand_min, 
		(const float)pso_rand_max
	);
	_pso_hw_real return_value = randGen.rand_num();
#else
	_pso_hw_real return_value = (_pso_hw_real)rand()/(_pso_hw_real)RAND_MAX;
#endif
	return return_value;
}

// ---------------------------------------------------
_pso_hw_real min_array(
	_pso_hw_real *array,//[_pso_n_S], 
	int * pos
){
	//[bestfitness(k), p] = min(f_ind);
#pragma HLS inline
	_pso_hw_real min = array[0];
	int pos_temp = 0;
	for (unsigned int i = 1; i < n_S; ++i) {
        if(array[i] < min) {
			min = array[i];
			pos_temp = i;
        }
	}
	*pos = pos_temp;
	return min;
}

// Init Functions
// ---------------------------------------------------
void initializeConstrains(
	volatile _pso_hw_real *u_curr
){
	// Initialize constrains based on current state
#pragma HLS inline
	initializeConstrains_loop: for (unsigned int i = 0; i < n_U; ++i){
		_pso_hw_real u_curr_tmp = u_curr[i];

        _pso_hw_real x_max_first_temp = u_curr_tmp + du_max[i]; //controled_system->acc_max[i];
		_pso_hw_real x_min_first_temp = u_curr_tmp + du_min[i];//controled_system->acc_min[i];        

		x_max_first[i] = (x_max_first_temp > u_max[i]) ? u_max[i]: x_max_first_temp;
		x_min_first[i] = (x_min_first_temp < u_min[i])? u_min[i]: x_min_first_temp;	
	}
}
// ---------------------------------------------------
void calculate_du_min(
){
#pragma HLS inline
	calculate_du_min_loop: for (unsigned i = 0; i < n_U; i++)
	{
		du_min[i] = -du_max[i];
	}
	
}
// ---------------------------------------------------
void initializeBestLocalFitness(void){
#pragma HLS inline
	// Initialize best local fitness 
	memset_loop<_pso_hw_real>(f_ind, (const _pso_hw_real)H_MAX, n_S);
	// for (unsigned int i = 0; i < n_S; ++i) { 
	// 	f_ind[i] = H_MAX;
	// }
}
// ---------------------------------------------------
int detectGlobalMinimum(
	int iter
){
	int pos;

	// Global Minimum detection
    bestfitness[iter] = min_array(f_ind, &pos);
    
	// // Global minimum position
    // for (unsigned int i = 0; i < Nu*n_U; ++i) {
    //     global_min[i] = _y[pos][i];
    // }
    return pos;

}
// ---------------------------------------------------
void equalizeParticles(
//	_pso_hw_real _x[][_pso_Nu*_pso_n_U],
	// _pso_hw_real **_x,
	_pso_hw_real iteration
){
//#pragma HLS inline
/*
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation
*/	
	_pso_hw_real x_tmp;
	if (iteration == 0){	
    for (unsigned int i = 0; i < n_S; ++i) {
// #pragma UNROLL
		for (unsigned int k = 1; k < Nu; ++k) {
	 	   	int idx1 = k*n_U;
	 	   	int idx2 = 0;
			for (unsigned int j = 0; j < n_U; ++j) {
				x_tmp = x[i*part_S + idx2];
                x[i*part_S + idx1] = x_tmp;
				idx1++;
				idx2++;
            }
        }
    }
	}
}
// ---------------------------------------------------
_pso_hw_real verifyControlConstrains(
	_pso_hw_real value, 
	int pos
){
#pragma HLS inline
	_pso_hw_real return_value;

	return_value = (value > u_max[pos]) ? u_max[pos] : value;
	return_value = (value < u_min[pos]) ? u_min[pos] : value;

	return return_value;
}
#ifdef PSO_CANON
// ---------------------------------------------------
void detectInvalidParticles(
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
){
// INITIALIZATION OF PARTICLES WITH Delta u CONTRAINS FOR MPC
//#pragma HLS inline
#pragma HLS ALLOCATION instances=rand_real limit=1 function
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation

	_pso_hw_real x_ant[n_U];

	initializeS_du_loop_top:for (unsigned int i = 0; i < n_S; ++i) {
		memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant, (_pso_hw_real *)u_curr);
		for (unsigned int j = 0; j < Nu; ++j) {
			for (unsigned int k = 0; k < n_U; ++k) {
#pragma HLS pipeline II=11 rewind
				int idx = j*n_U + k;
				_pso_hw_real rand_tmp = rand_real();
				_pso_hw_real du_tmp = rand_tmp*(_pso_hw_real)2.0*du_max[k] - du_max[k];
		        _pso_hw_real x_tmp = x_ant[k] + du_tmp; //random->read(); 
		        // _pso_hw_real x_tmp = ( x_ant[k] + (-du_max[k]) ) + ((_pso_hw_real)2.0*du_max[k]) * rand_tmp; //random->read(); 
				x_tmp = verifyControlConstrains(x_tmp, k);
                x[i*part_S + idx] = x_tmp;
				y[i*part_S + idx] = x_tmp; // uss_local[k];
		        v[i*part_S + idx] = init_v;
                valid_particle[i] = 1;
				x_ant[k] = x_tmp;
				idx++;
		    }
		}
	}
}

#ifdef PSO_CANON
// ---------------------------------------------------
void initializeParticles(
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
				_pso_hw_real rand_tmp = rand_real();
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
void initializeLastBestKPSO(
	_pso_hw_real volatile *last_best,
	// _pso_hw_real _x[][_pso_Nu*_pso_n_U]
	// _pso_hw_real **_x,
	uint8_t index
){
// #pragma HLS inline
	// Uses KPSO (best last position)

	for (unsigned int i = 0; i < n_U; i++){
		_pso_hw_real last_best_tmp = last_best[i];
		_pso_hw_real x_local;
		x_local = (last_best_tmp > x_max_first[i]) ? x_max_first[i] : last_best_tmp;
		x_local = (last_best_tmp < x_min_first[i]) ? x_max_first[i] : last_best_tmp;
		x[index*part_S + i] = x_local;
	}
	
	for (unsigned int i = n_U; i < (Nu*n_U); i++) {
		x[index*part_S + i] = last_best[i];			
	}
}
// ---------------------------------------------------
void initializeStableZero(
//	_pso_hw_real *uss, 
	_pso_hw_real *u_curr, 	
	// _pso_hw_real _x[][_pso_Nu*_pso_n_U]
	// _pso_hw_real **_x,
	uint8_t index
){
// #pragma HLS inline
// Starts one particle with all Zeros for 
// stable response after equilibrium is reached
    int idx;
	_pso_hw_real x_ant[_pso_n_U];
	memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_n_U>(x_ant, (_pso_hw_real*)u_curr);
   	for (unsigned int k = 0; k < Nu; ++k){
        for (unsigned int i = 0; i < n_U; ++i) {
        	idx = k*n_U + i;
		    _pso_hw_real x_tmp = uss_local[i];
			// _pso_hw_real comp_tmp = (k == 0) ? x_tmp - u_curr[i] : x_tmp - x_ant[i];
			_pso_hw_real comp_tmp = x_tmp - x_ant[i] ;
			
			x_tmp = (comp_tmp > du_max[i]) ? du_max[i] : x_tmp;
			x_tmp = (comp_tmp < du_min[i]) ? du_min[i] : x_tmp;

            x[index*part_S + idx] = verifyControlConstrains(x_tmp, i);
			x_ant[i] = x_tmp;
			idx++;
        }
   	}
}
// ---------------------------------------------------
void evaluateFitnessAndDetectLocalBest(
	// _pso_hw_real _x[][_pso_Nu*_pso_n_U], 
	// _pso_hw_real _y[][_pso_Nu*_pso_n_U],
	// _pso_hw_real **_x,
	// _pso_hw_real **_y,

	_pso_hw_real *x_curr_local,//[_pso_Nx],

	_pso_hw_real *xref,//[_pso_Nx*_pso_Nu], 
	_pso_hw_real *uref,//[_pso_n_U],
	_pso_hw_real *xss//[_pso_Nx]
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
		current_system.nmpc_cost_function_topflow(x_curr_local, &x[i*part_S], xref, &fx[i]);
		// fx[i] = rand_real();
		if (fx[i] < f_ind[i]) {
			memcpy_loop_rolled<_pso_hw_real, _pso_hw_real, _pso_Nu*_pso_n_U>(&y[i*part_S], &x[i*part_S]);
			f_ind[i] = fx[i];           
		}
        //}
    }
    // std::cout << std::endl;

}
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
){    

    // Update Particles
    for (unsigned int i = 0; i < n_S; ++i) {
		for (unsigned int k = 0; k < Nu; ++k) {
	 	   	for (unsigned int j = 0; j < n_U; ++j) {
    	 		int idx = k*Nu+j;
	
				_pso_hw_real v_tmp = v[i*part_S + idx];
				_pso_hw_real x_tmp = x[i*part_S + idx];
				
				_pso_hw_real r1 = rand_real(); //random->read();
				_pso_hw_real r2 = rand_real(); //random->read();

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
void updateParticlesWithDuConstrains(
//	_pso_hw_real global_min[],

	// _pso_hw_real *_x[][_pso_Nu*_pso_n_U], 
	// _pso_hw_real *_y[][_pso_Nu*_pso_n_U], 
	// _pso_hw_real *_v[][_pso_Nu*_pso_n_U]
	// _pso_hw_real **_x,
	// _pso_hw_real **_y,
	// _pso_hw_real **_v
){
#pragma HLS ALLOCATION type=function instances=rand_real limit=1 

	// Update particles
	_pso_hw_real x_ant[n_U];
    for (unsigned int i = 0; i < n_S; ++i) {
// #pragma UNROLL
#pragma HLS ALLOCATION type=operation instances=hmul 	  limit=1 
#pragma HLS ALLOCATION type=operation instances=hadd 	  limit=1 
#pragma HLS ALLOCATION type=operation instances=hsub 	  limit=1 
#pragma HLS ALLOCATION type=function  instances=rand_real limit=1 
		memset_loop<_pso_hw_real>(x_ant, (const _pso_hw_real)0.0, n_U);
        //if(valid_particle[i] == 1){
    	 for (unsigned int j = 0; j < Nu; ++j) {
    	 	for (unsigned int k = 0; k < n_U; ++k){
#pragma HLS pipeline II=11 rewind
    	 		// int idx = k*Nu+j;
				int idx = j*n_U+k;
                _pso_hw_real r1 = rand_real(); //random->read();
                _pso_hw_real r2 = rand_real(); //random->read();

	            _pso_hw_real v_tmp = v[i*part_S + idx];
	            _pso_hw_real x_tmp = x[i*part_S + idx];
	            // v = w*v + c1*r1*(y-x) + c2*r2*(global_min - x)
				_pso_hw_real v_new = w*v_tmp + c1*r1*(y[i*part_S + idx]-x_tmp) + c2*r2*(global_min[idx] - x_tmp);
	            v_new = (v_new > max_v) ? max_v : v_new;
	            v_new = (v_new < min_v) ? min_v : v_new;
				v[i*part_S + idx] = v_new;
	            
				_pso_hw_real x_new = x_tmp + v_new;
				_pso_hw_real x_tmp2 = 	(j==0) ? x_new : x_new - x_ant[k]; //x[i][idx]-x[i][idx-1];
				_pso_hw_real cmp_max = 	(j==0) ? x_max_first[k] : du_max[k];
				_pso_hw_real cmp_min = 	(j==0) ? x_min_first[k] : du_min[k];
				x_new = (x_tmp2 > cmp_max) ? x_ant[k] + cmp_max : x_new;
				x_new = (x_tmp2 < cmp_min) ? x_ant[k] + cmp_min : x_new;
                x[i*part_S + idx] = verifyControlConstrains(x_new, k);

				x_ant[k] = x_new;
	            
            } // END FOR k
        } // END FOR j
       //} // END IF
    } // END FOR i

}

	
};

#endif