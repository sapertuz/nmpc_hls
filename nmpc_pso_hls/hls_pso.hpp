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
//     #define _Nuu    _N
//     #define _n_U     1
// #else
//     #define _Nn     _N
//     #define _Nuu    _Nu
//     #define _n_U     _n_U
// #endif 

template<
    class _hw_real,
    class _randCore_t,
    class _system_t,
	unsigned _n_S,
	unsigned _maxiter,
    unsigned _Nh,
    unsigned _Nx,
    unsigned _n_U,
    unsigned _Nu
>class PSO{
protected:
	// Pseudo random generator variables
	

	const unsigned int n_S = _n_S; 	// number of particles
	//System * controled_system;

	// const int kpso = _KPSO;
	const int stable_zero;
	const uint8_t particle_last_best = 0;
	const uint8_t particle_stable_zero = 1;

	// PSO Parameters
	const unsigned int	maxiter = _maxiter;
	const _hw_real 	max_v;
	const _hw_real 	min_v;
	const _hw_real 	w0; // initial weight
	const _hw_real	wf; // final weight
	const _hw_real	slope;
	const _hw_real	c1; // cognitive coefficient
	const _hw_real	c2; // social coefficient
	
	_hw_real w;
	const _hw_real slope_init;
	
	_hw_real threshold;
	int stop_criteria;

	const _hw_real *u_min;//[_n_U];
	const _hw_real *u_max;//[_n_U];
	const _hw_real *du_max;//[_n_U];
	_hw_real du_min[_n_U];//[_n_U];
	const _hw_real *uss_local;//[_n_U];
	_hw_real x_max_first[_n_U];
	_hw_real x_min_first[_n_U];

	const _hw_real init_v; // initial velocity

	// Particles
	// _hw_real x[_n_S][_Nu*_n_U];
	// _hw_real y[_n_S][_Nu*_n_U];
	// _hw_real v[_n_S][_Nu*_n_U];
	int valid_particle[_n_S];
	int number_of_active_particles;

	_hw_real f_ind[_n_S];
	_hw_real fx[_n_S];
	_hw_real bestfitness[_maxiter];
	_hw_real global_min[_Nu*_n_U];
	//_hw_real previous_control[_Parametrization];

	// Controlled System Configuration
	const unsigned int Nh = _Nh;
	const unsigned int Nx = _Nx;
	const unsigned int n_U = _n_U;
	const unsigned int Nu = _Nu;
	// _hw_real u_from_parameters[_N*_n_U];

	//ExportData * exporter;

/*	
	// System
	const unsigned short controlled_state[_Nx];
	const _hw_real state_upper_limits[_Nx];
	const _hw_real state_lower_limits[_Nx];
	const _hw_real Q[_Nx];
	const _hw_real Qf[_Nx];
	const _hw_real R[_n_U];
	const _hw_real Ts;
	
	typedef System<_hw_real,_Nh, _Nx, _n_U, _Nu> _system_t; 
	

	// Random gen
	const _hw_real rand_min = -1.0;
	const _hw_real rand_max = 1.0;
	const int rand_seed[_n_S];

	typedef pseudoRand_gen<_hw_real, _n_S> _randCore_t;
*/
	_system_t *current_system;
	_randCore_t *randGen;
public:
// Init Execute
constexpr PSO(
	// PSO Configuration
	const int 	 	_stable_zero,
	const _hw_real  _max_v,
	const _hw_real  _w0,
	const _hw_real 	_wf,
	const _hw_real 	_slope,
	const _hw_real 	_c1,
	const _hw_real 	_c2,
	const _hw_real  _slope_init,
	const _hw_real	*_u_min,
	const _hw_real	*_u_max,
	const _hw_real	*_du_max,
	const _hw_real 	*_uss
/*	
	,
	const unsigned short _controlled_state[_Nx],
	const _hw_real _state_upper_limits[_Nx],
	const _hw_real _state_lower_limits[_Nx],
	const _hw_real _Q[_Nx],
	const _hw_real _Qf[_Nx],
	const _hw_real _R[_n_U],
	const _hw_real _Ts,
	const int _rand_seed[_n_S]
*/
	,
	_system_t *_current_system,
	_randCore_t *_randGen
	) : max_v(_max_v), min_v(-_max_v), w0(_w0), wf(_wf), 
		slope(_slope), c1(_c1), c2(_c2), slope_init((_wf-_w0)/_maxiter),
		stable_zero(_stable_zero), init_v(_max_v*0.01),
		u_min(_u_min), u_max(_u_max), du_max(_du_max), uss_local(_uss),
/*		
		Ts(_Ts), controlled_state(_controlled_state), 
        state_upper_limits(_state_upper_limits), state_lower_limits(_state_lower_limits),
        Q(_Q), Qf(_Qf), R(_R), 
*/
		randGen(_randGen), current_system(_current_system)
{
}

// ---------------------------------------------------
int execute(
	volatile _hw_real x_curr[_Nx], 
	volatile _hw_real u_curr[_n_U], 
	int iteration, 
	volatile _hw_real last_best[_n_U*_Nu], 
	volatile _hw_real xref[_Nx*_Nh],
	volatile _hw_real uref[_n_U], 
	volatile _hw_real xss[_Nx],
	//volatile _hw_real uss[_n_U], 
	
	_hw_real new_best[_Nu*_n_U],
	_hw_real * J
){

#pragma HLS ALLOCATION operation instances=hadd limit=2
#pragma HLS ALLOCATION operation instances=hsub limit=2
#pragma HLS ALLOCATION operation instances=hmul limit=2

#pragma HLS bind_storage variable=x_max_first type=RAM_2P impl=LUTRAM
#pragma HLS bind_storage variable=x_min_first type=RAM_2P impl=LUTRAM
    // Update sensor read
    //static System * controled_system = new ModelState();
    //controled_system->setState(x_curr);

	// Particles Variables
#ifdef __SYNTHESIS__
	_hw_real x[n_S][Nu*n_U];
#pragma HLS bind_storage variable=x type=RAM_T2P impl=BRAM
#pragma HLS stream variable=x type=shared
	_hw_real y[n_S][Nu*n_U];
#pragma HLS bind_storage variable=y type=RAM_T2P impl=BRAM
#pragma HLS stream variable=y type=shared
	_hw_real v[n_S][Nu*n_U];
#pragma HLS bind_storage variable=v type=RAM_T2P impl=BRAM
#pragma HLS stream variable=v type=shared
#else
	//_hw_real bestfitness[_maxiter];
	//_hw_real global_min[Nu*n_U];
	_hw_real **x;
	_hw_real **y;
	_hw_real **v;
	if((x = alloc_matrix(n_S, Nu*n_U)) == NULL) {return -1;}
	if((y = alloc_matrix(n_S, Nu*n_U)) == NULL) {return -1;}
	if((v = alloc_matrix(n_S, Nu*n_U)) == NULL) {return -1;}

#endif

    number_of_active_particles = n_S;    
	int best_pos;

#ifdef __SYNTHESIS__
	_hw_real u_curr_local[n_U];
#pragma HLS bind_storage variable=u_curr_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=u_curr_local type=shared

	_hw_real x_curr_local[Nx];
#pragma HLS bind_storage variable=x_curr_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=x_curr_local type=shared

	_hw_real xref_local[Nx*Nh];
#pragma HLS bind_storage variable=xref_local type=RAM_2P impl=BRAM
#pragma HLS stream variable=xref_local type=shared

	_hw_real uref_local[n_U];
#pragma HLS bind_storage variable=uref_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=uref_local type=shared

	_hw_real xss_local[Nx];
#pragma HLS bind_storage variable=xss_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=xss_local type=shared
#else
	_hw_real *u_curr_local;
	_hw_real *x_curr_local;
	_hw_real *xref_local;
	_hw_real *uref_local;
	_hw_real *xss_local;
	x_curr_local = (_hw_real *) malloc(Nx*sizeof(_hw_real));
	u_curr_local = (_hw_real *) malloc(n_U*sizeof(_hw_real));
	xref_local = (_hw_real *) malloc(Nx*Nh*sizeof(_hw_real));
	uref_local = (_hw_real *) malloc(n_U*sizeof(_hw_real));
	xss_local = (_hw_real *) malloc(Nx*sizeof(_hw_real));
#endif
	memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(u_curr_local, u_curr);
	memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(x_curr_local, xss);
	memcpy_loop_rolled<_hw_real, _hw_real, _Nx*_Nh>(xref_local, xref);
	memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(uref_local, uref);
	memcpy_loop_rolled<_hw_real, _hw_real, _Nx>(xss_local, xss);

/*
	_hw_real uss_local[n_U];
#pragma HLS bind_storage variable=uss_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=uss_local type=shared
	memcpy_loop_rolled<_hw_real, _hw_real, n_U>(uss_local, uss);
*/
	calculate_du_min();

	initializeConstrains(u_curr_local);
    initializeParticlesWithDuConstrains(u_curr_local, x, y, v);

//#if (n_U > 1)
//    if((iteration == 0) && (n_U > 1)) {
//	equalizeParticles(x, iteration);
//    }        
//#endif

//    if ((kpso == 1) && (iteration > 0)) {
    initializeLastBestKPSO(last_best, x, particle_last_best);
//    }
//    if (stable_zero == 1) {
	initializeStableZero(u_curr_local, x, particle_stable_zero);
//    } 

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
			x, 
			y, 
			x_curr_local, 
			xref_local, 
			uref_local, 
			xss_local
		);

        // Global Minimum detection
		best_pos = detectGlobalMinimum(k);
		memcpy_loop_rolled<_real, _real, _Nu*_n_U>(global_min, &y[best_pos][0]);

        updateParticlesWithDuConstrains(x, y, v);    
#ifdef DEBUG_PSO
		std::cout << "[" << k << "] ";
		std::cout << "\t"; print_formatted_float_array(&global_min[0], n_U, 2, 6);
		std::cout << "|";
		std::cout << "\t"; print_formatted_float_array(f_ind, n_S, 2, 6); 
		std::cout << "|";
		std::cout << std::endl;
#endif
        //detectInvalidParticles(k, best_pos);
/*
		// STOP CRITERIA
		if(stop_criteria) {
			if (k > stop_criteria) {
			    if ((bestfitness[k-stop_criteria]-bestfitness[k]) < threshold)
			        break;
			}
		}
        if(number_of_active_particles < 2){
            break;
        }
*/
		k++;
		w+= slope;
	}
    
	// Return best value
    for (unsigned int i = 0; i < Nu*n_U; ++i) {
        new_best[i] = global_min[i];
    }
    
    J[0] = bestfitness[_maxiter-1];

#ifndef __SINTHESYS__
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

// Misc Functions
// ---------------------------------------------------
_hw_real rand_real(unsigned core){
#ifndef __SYNTHESIS__
	_hw_real return_value = (_hw_real)rand()/(_hw_real)RAND_MAX;
#else
	_hw_real return_value = randGen->rand_num(core);
#endif
	return return_value;
}
_hw_real rand_real(){
#ifndef __SYNTHESIS__
	_hw_real return_value = (_hw_real)rand()/(_hw_real)RAND_MAX;
#else
	_hw_real return_value = randGen->rand_num(0);
#endif
	return return_value;
}

// ---------------------------------------------------
_hw_real min_array(
	_hw_real *array,//[_n_S], 
	int * pos
){
	//[bestfitness(k), p] = min(f_ind);
#pragma HLS inline
	_hw_real min = array[0];
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
	volatile _hw_real *u_curr
){
	// Initialize constrains based on current state
#pragma HLS inline
	for (unsigned int i = 0; i < n_U; ++i){
		_hw_real u_curr_tmp = u_curr[i];

        _hw_real x_max_first_temp = u_curr_tmp + du_max[i]; //controled_system->acc_max[i];
		x_max_first[i] = (x_max_first_temp > u_max[i]) ? u_max[i]: x_max_first_temp;

		_hw_real x_min_first_temp = u_curr_tmp + (-du_max[i]);//controled_system->acc_min[i];        
		x_min_first[i] = (x_min_first_temp < u_min[i])? u_min[i]: x_min_first_temp;	
	}
}
// ---------------------------------------------------
void calculate_du_min(
){
	for (unsigned i = 0; i < n_U; i++)
	{
		du_min[i] = -du_max[i];
	}
	
}
// ---------------------------------------------------
void initializeBestLocalFitness(void){
#pragma HLS inline
	// Initialize best local fitness 
	memset_loop<_hw_real>(f_ind, (const _hw_real)H_MAX, n_S);
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
//	_hw_real _x[][_Nu*_n_U],
	_hw_real **_x,
	_hw_real iteration
){
//#pragma HLS inline
/*
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation
*/	
	_real x_tmp;
	if (iteration == 0){	
    for (unsigned int i = 0; i < n_S; ++i) {
#pragma UNROLL
		for (unsigned int k = 1; k < Nu; ++k) {
	 	   	int idx1 = k*n_U;
	 	   	int idx2 = 0;
			for (unsigned int j = 0; j < n_U; ++j) {
				x_tmp = _x[i][idx2];
                _x[i][idx1] = x_tmp;
				idx1++;
				idx2++;
            }
        }
    }
	}
}
// ---------------------------------------------------
_hw_real verifyControlConstrains(
	_hw_real value, 
	int pos
){
#pragma HLS inline
	_hw_real return_value;

	return_value = (value > u_max[pos]) ? u_max[pos] : value;
	return_value = (value < u_min[pos]) ? u_min[pos] : value;

	return return_value;
}

// ---------------------------------------------------
void detectInvalidParticles(
	int iter, 
	int best_pos,
	
//	_hw_real _x[][_Nu*_n_U]
	_hw_real *_x
){
	for (unsigned int i = 0; i < n_S; ++i) {
        if(n_S != best_pos) {
            if((std::abs(fx[i]-bestfitness[iter])/bestfitness[iter]) < 0.0001){
//                int count = 0;  
                int count = Nh*n_U;
                for (unsigned int j = 0; j < Nh*n_U; ++j) {
                    if((std::abs(_x[i][j]-global_min[j])/global_min[j]) < 0.0001){
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

// Workflow Functions
// ---------------------------------------------------
void initializeParticlesWithDuConstrains(
	volatile _hw_real *u_curr, 
	//volatile _hw_real *uss,

	// _hw_real _x[][_Nu*_n_U], 
	// _hw_real _y[][_Nu*_n_U], 
	// _hw_real _v[][_Nu*_n_U]
	_hw_real ** _x,
	_hw_real ** _y,
	_hw_real ** _v
){
// INITIALIZATION OF PARTICLES WITH Delta u CONTRAINS FOR MPC
#pragma HLS inline
#pragma HLS ALLOCATION instances=rand_real limit=1 function
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation

	_hw_real x_ant[n_U];
	memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(x_ant, u_curr);

	for (unsigned int i = 0; i < n_S; ++i) {
#pragma UNROLL
		for (unsigned int j = 0; j < Nu; ++j) {
			int idx = j*n_U;
			for (unsigned int k = 0; k < n_U; ++k) {
#pragma HLS pipeline II=11 rewind
				_hw_real rand_tmp = rand_real();
		        _hw_real x_tmp = (x_ant[k] + (-du_max[k])) + ((_hw_real)2.0*du_max[k]) * rand_tmp; //random->read(); 
				_x[i][idx] = verifyControlConstrains(x_tmp, k);
                _y[i][idx] = uss_local[k];
		        _v[i][idx] = init_v;
                valid_particle[i] = 1;
				x_ant[k] = x_tmp;
				idx++;
		    }
		}
	}
}

// ---------------------------------------------------
void initializeParticles(
    // _hw_real _x[][_Nu*_n_U], 
    // _hw_real _y[][_Nu*_n_U], 
    // _hw_real _v[][_Nu*_n_U]
	_hw_real **_x, 
    _hw_real **_y, 
    _hw_real **_v
) {
// GENERAL PSO INITIALIZATION OF PARTICLES
	_hw_real d_max = 0.2;
    _hw_real d2_max = 0.4;
	_hw_real x_ant[n_U];
	for (unsigned int i = 0; i < n_S; ++i) {
#pragma UNROLL
		for (unsigned int j = 0; j < Nu; ++j) {
			int idx = j*n_U;
			for (unsigned int k = 0; k < n_U; ++k) {
				_hw_real rand_tmp = rand_real();
				_hw_real x_tmp = (j == 0)? 
					u_min[k] + (u_max[k]-u_min[k]) * rand_tmp: //random->read();
					x_ant[k] - d_max + (d2_max) * rand_tmp; //random->read();
				_x[i][idx] = x_tmp; 
				_y[i][idx] = H_MAX;
				_v[i][idx] = init_v;  
				x_ant[k] = x_tmp;
			}
        }
        valid_particle[i] = 1;
    }
}

// ---------------------------------------------------
void initializeLastBestKPSO(
	_hw_real volatile *last_best,
	// _hw_real _x[][_Nu*_n_U]
	_hw_real **_x,
	uint8_t index
){
#pragma HLS inline
	// Uses KPSO (best last position)

	for (unsigned int i = 0; i < n_U; i++){
		_hw_real last_best_tmp = last_best[i];
		_hw_real x_local;
		x_local = (last_best_tmp > x_max_first[i]) ? x_max_first[i] : last_best_tmp;
		x_local = (last_best_tmp < x_min_first[i]) ? x_max_first[i] : last_best_tmp;
		_x[index][i] = x_local;
	}
	
	for (unsigned int i = n_U; i < (Nu*n_U); i++) {
		_x[index][i] = last_best[i];			
	}
}
// ---------------------------------------------------
void initializeStableZero(
//	_hw_real *uss, 
	_hw_real *u_curr, 	
	// _hw_real _x[][_Nu*_n_U]
	_hw_real **_x,
	uint8_t index
){
#pragma HLS inline
// Starts one particle with all Zeros for 
// stable response after equilibrium is reached
    int idx;
	_hw_real x_ant[n_U];
	memset_loop<_hw_real>(x_ant, (const _hw_real)0.0, n_U);
   	for (unsigned int k = 0; k < Nu; ++k){
        idx = k*n_U;
		for (unsigned int i = 0; i < n_U; ++i) {
            _hw_real x_tmp = uss_local[i];
			_hw_real comp_tmp = (k == 0) ? x_tmp - u_curr[i] : x_tmp - x_ant[i];

			x_tmp = (comp_tmp > du_max[i]) ? du_max[i] : x_tmp;
			x_tmp = (comp_tmp < du_min[i]) ? du_min[i] : x_tmp;

            _x[index][idx] = verifyControlConstrains(x_tmp, i); 
			x_ant[i] = x_tmp;
			idx++;
        }
   	}
}
// ---------------------------------------------------
void evaluateFitnessAndDetectLocalBest(
	// _hw_real _x[][_Nu*_n_U], 
	// _hw_real _y[][_Nu*_n_U],
	_hw_real **_x,
	_hw_real **_y,

	_hw_real *x_curr_local,//[_Nx],

	_hw_real *xref,//[_Nx*_Nu], 
	_hw_real *uref,//[_n_U],
	_hw_real *xss//[_Nx]
){

// #pragma HLS interface ap_bus  port=_x
// #pragma HLS interface ap_bus  port=xref
// #pragma HLS interface ap_fifo port=uref
// #pragma HLS interface ap_fifo port=xss

/*
	typedef System<_real, _Nh, _Nx, _n_U, _Nu> T_system;
    T_system current_system(
        u_max, 
        u_min, 
        du_max,
        controlled_state,
        state_upper_limits, 
        state_lower_limits, 
        Q, 
        Qf, 
        R, 
        uss, 
        Ts);
*/

	// Evaluates fitness and local detection
    loop_1: for(unsigned int i = 0; i < n_S; i++) {
        // std::cout << std::endl;
#pragma UNROLL
#pragma HLS unroll factor=4 region skip_exit_check
        //if(valid_particle[i] == 1){
            // std::cout << x[i][0] << " " << x[i][1] << " " << x[i][2] << " " << x[i][3];
            // fx[i] = rand_real(); 
			current_system->nmpc_cost_function(x_curr_local, &_x[i][0], xref, &fx[i]);
            if (fx[i] < f_ind[i]) {
				memcpy_loop_rolled<_hw_real, _hw_real, _Nu*_n_U>(&_y[i][0], &_x[i][0]);
                // loop_2: for (unsigned int j = 0; j < _Nu*n_U; ++j) {
                // 	_y[i][j] = _x[i][j];
                // }
                f_ind[i] = fx[i];           
            }
        //}
    }
    // std::cout << std::endl;

}
// ---------------------------------------------------
void updateParticles(
//	_hw_real global_min[],

    // _hw_real _x[][_Nu*_n_U],
    // _hw_real _y[][_Nu*_n_U],
    // _hw_real _v[][_Nu*_n_U]
	_hw_real **_x,
	_hw_real **_y,
	_hw_real **_v
){
    _hw_real r1, r2;
    
    // Update Particles
    for (unsigned int i = 0; i < n_S; ++i) {
#pragma UNROLL
		for (unsigned int k = 0; i < n_U; ++i) {
	 	   for (unsigned int j = 0; j < Nu; ++j) {
    	 		int idx = k*Nu+j;
                r1 = rand_real();
                r2 = rand_real();
	
				_hw_real v_tmp = _v[i][idx];
				_hw_real x_tmp = _x[i][idx];
				// v = w*v + c1*r1*(y-x) + c2*r2*(global_min - x)
				_hw_real v_new = w*v_tmp + c1*r1*(_y[i][idx]-x_tmp) + c2*r2*(global_min[idx] - x_tmp);

                v_new = (v_new > max_v) ? max_v : v_new;
				v_new = (v_new < min_v) ? min_v : v_new;
				_v[i][idx] = v_new;

				_hw_real x_new = x_tmp + v_new;
				_x[i][idx] = verifyControlConstrains(x_new, k);
            } // END FOR j
        } // End FOR k
    } // END FOR i
}
// ---------------------------------------------------
void updateParticlesWithDuConstrains(
//	_hw_real global_min[],

	// _hw_real *_x[][_Nu*_n_U], 
	// _hw_real *_y[][_Nu*_n_U], 
	// _hw_real *_v[][_Nu*_n_U]
	_hw_real **_x,
	_hw_real **_y,
	_hw_real **_v
){
#pragma HLS interface ap_bus  port=x
#pragma HLS interface ap_bus  port=y
#pragma HLS interface ap_bus  port=v

#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation

#pragma HLS ALLOCATION instances=rand_real limit=1 function

    _hw_real r1, r2;
	_hw_real x_ant[n_U];
	_hw_real v_tmp = 0.0;
	_hw_real x_tmp = 0.0;
	_hw_real v_new = 0.0;
	_hw_real x_new = 0.0;
	//memset(x_ant, (const _hw_real)0.0, n_U*sizeof(_hw_real));
    
	// First Moment of each Particle
	for (unsigned int i = 0; i < n_S; ++i) {
#pragma UNROLL
		for (unsigned int k = 0; k < n_U; ++k){
			int idx = k*Nu;
			r1 = rand_real(); //random->read();
			r2 = rand_real(); //random->read();

			v_tmp = _v[i][idx];
			x_tmp = _x[i][idx];
			// v = w*v + c1*r1*(y-x) + c2*r2*(global_min - x)
			v_new = w*v_tmp + c1*r1*(_y[i][idx]-x_tmp) + c2*r2*(global_min[idx] - x_tmp);
			v_new = (v_new > max_v) ? max_v : v_new;
			v_new = (v_new < min_v) ? min_v : v_new;
			_v[i][idx] = v_new;
			
			x_new = x_tmp + v_new;
			x_new = (x_new > x_max_first[k]) ? x_max_first[k] : x_new;
			x_new = (x_new < x_min_first[k]) ? x_min_first[k] : x_new;
			_x[i][idx] = verifyControlConstrains(x_new, k);

			x_ant[k] = x_new;
			
			
		} // END FOR k
    } // END FOR i

	// Update rest of particles
    for (unsigned int i = 0; i < n_S; ++i) {
#pragma UNROLL
        //if(valid_particle[i] == 1){
    	 for (unsigned int j = 1; j < Nu; ++j) {
    	 	for (unsigned int k = 0; k < n_U; ++k){
    	 		int idx = k*Nu+j;
                r1 = rand_real(); //random->read();
                r2 = rand_real(); //random->read();

	            _hw_real v_tmp = _v[i][idx];
	            _hw_real x_tmp = _x[i][idx];
	            // v = w*v + c1*r1*(y-x) + c2*r2*(global_min - x)
				_hw_real v_new = w*v_tmp + c1*r1*(_y[i][idx]-x_tmp) + c2*r2*(global_min[idx] - x_tmp);
	            v_new = (v_new > max_v) ? max_v : v_new;
	            v_new = (v_new < min_v) ? min_v : v_new;
				_v[i][idx] = v_new;
	            
				_hw_real x_new = x_tmp + v_new;
				_hw_real x_tmp2 = x_tmp - x_ant[k]; //x[i][idx]-x[i][idx-1];
				x_new = (x_tmp2 > du_max[k]) ? x_ant[k] + du_max[k] : x_new;
				x_new = (x_tmp2 < du_min[k]) ? x_ant[k] - du_max[k] : x_new;
                _x[i][idx] = verifyControlConstrains(x_new, k);

				x_ant[k] = x_new;
	            
            } // END FOR k
        } // END FOR j
       //} // END IF
    } // END FOR i
}
};

#endif