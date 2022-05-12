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
	unsigned _n_S,
	unsigned _maxiter,
    unsigned _Nh,
    unsigned _Nx,
    unsigned _n_U,
    unsigned _Nu
>class PSO{
protected:
	// Pseudo random generator variables
	

	const int S; 	// number of particles
	//System * controled_system;

	// const int kpso = _KPSO;
	const int stable_zero;

	// PSO Parameters
	const int 		maxiter;
	const _hw_real 	max_v;
	const _hw_real 	w0; // initial weight
	const _hw_real	wf; // final weight
	const _hw_real	slope;
	const _hw_real	c1; // cognitive coefficient
	const _hw_real	c2; // social coefficient
	
	_hw_real w;
	const _hw_real slope_init;
	
	_hw_real threshold;
	int stop_criteria;

	const _hw_real u_min[_n_U];
	const _hw_real u_max[_n_U];
	const _hw_real du_max[_n_U];
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
	int N;
	int Nu;
	int Nc;

	// _hw_real u_from_parameters[_N*_n_U];

	//ExportData * exporter;
	
	// System
	const unsigned short controlled_state[_Nx];
	const _hw_real state_upper_limits[_Nx];
	const _hw_real state_lower_limits[_Nx];
	const _hw_real Q[_Nx];
	const _hw_real Qf[_Nx];
	const _hw_real R[_n_U];
	const _hw_real uss[_n_U];
	const _hw_real Ts;
	
	// typedef System<_hw_real,_Nh, _Nx, _n_U, _Nu> _system_t; 
	// _system_t current_system;
	typedef pseudoRand_gen<_hw_real> _randCore_t;
	_randCore_t randGen(seed, min, max);

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
	const _hw_real	_u_min[_n_U],
	const _hw_real	_u_max[_n_U],
	const _hw_real	_du_max[_n_U],
	const unsigned short _controlled_state[_Nx],
	const _hw_real _state_upper_limits[_Nx],
	const _hw_real _state_lower_limits[_Nx],
	const _hw_real _Q[_Nx],
	const _hw_real _Qf[_Nx],
	const _hw_real _R[_n_U],
	const _hw_real _uss[_n_U],
	const _hw_real _Ts
	) : max_v(_max_v), w0(_w0), wf(_wf), 
		slope(_slope), c1(_c1), c2(_c2), slope_init((_wf-_w0)/_maxiter),
		stable_zero(_stable_zero), init_v(_max_v*0.01),
		u_min(_u_min), u_max(_u_max), du_max(_du_max),
		Ts(_Ts), uss(_uss), controlled_state(_controlled_state), 
        state_upper_limits(_state_upper_limits), state_lower_limits(_state_lower_limits),
        Q(_Q), Qf(_Qf), R(_R)
{
}

// ---------------------------------------------------
int execute(
	volatile _hw_real x_curr[_Nx], 
	volatile _hw_real u_curr[_n_U], 
	int iteration, 
	volatile _hw_real last_best[_Nu*_n_U], 
	volatile _hw_real xref[_Nu*_Nx],
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
	_hw_real x[_n_S][_Nu*_n_U];
#pragma HLS bind_storage variable=x type=RAM_T2P impl=BRAM
#pragma HLS stream variable=x type=shared
	_hw_real y[_n_S][_Nu*_n_U];
#pragma HLS bind_storage variable=y type=RAM_T2P impl=BRAM
#pragma HLS stream variable=y type=shared
	_hw_real v[_n_S][_Nu*_n_U];
#pragma HLS bind_storage variable=v type=RAM_T2P impl=BRAM
#pragma HLS stream variable=v type=shared
	//_hw_real bestfitness[_maxiter];
	//_hw_real global_min[_Nu*_n_U];
    number_of_active_particles = _n_S;
    
	int best_pos;

	_hw_real u_curr_local[_n_U];
#pragma HLS bind_storage variable=u_curr_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=u_curr_local type=shared
	memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(u_curr_local, u_curr);

	_hw_real uss_local[_n_U];
#pragma HLS bind_storage variable=uss_local type=RAM_2P impl=LUTRAM
#pragma HLS stream variable=uss_local type=shared
	memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(uss_local, uss);

	initializeConstrains(u_curr_local);
    initializeParticlesWithDuConstrains(u_curr_local, uss_local, x, y, v);

//#if (_n_U > 1)
//    if((iteration == 0) && (_n_U > 1)) {
	equalizeParticles(x, iteration);
//    }        
//#endif

//    if ((kpso == 1) && (iteration > 0)) {
    initializeLastBestKPSO(last_best, x);
//    }
    if (stable_zero == 1) {
        initializeStableZero(uss_local, u_curr_local, 1, x);
    } 

    initializeBestLocalFitness();
	// ITERATIVE PROCESS
	int k = 0;  // index of iteration

	while (k < _maxiter) {
        evaluateFitnessAndDetectLocalBest(x, y, xref, uref, xss, uss_local);     
        // Global Minimum detection
		best_pos = detectGlobalMinimum(k);
		memcpy_loop_rolled<_real, _real, _Nu*_n_U>(global_min, y[best_pos]);

        updateParticlesWithDuConstrains(x, y, v);    

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
		k = k + 1;
		w = w + slope;
	}
    
	// Return best value
    for (int i = 0; i < _Nu*_n_U; ++i) {
        new_best[i] = global_min[i];
    }
    
    *J = bestfitness[_maxiter-1];

	return 0;
}

// Misc Functions
// ---------------------------------------------------
_hw_real rand_real(){
#ifndef __SYNTHESIS__
	return (_hw_real)rand()/(_hw_real)RAND_MAX;
#else
	_hw_real return_value = randGen.rand_num();
    return return_value;
#endif
}

// ---------------------------------------------------
_hw_real min_array(
	_hw_real array[_n_S], 
	int * pos
){
	//[bestfitness(k), p] = min(f_ind);
#pragma HLS inline
	_hw_real min = array[0];
	int pos_temp = 0;
	for (int i = 1; i < _n_S; ++i) {
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
	for (int i = 0; i < _n_U; ++i){
		_hw_real u_curr_local = u_curr[i];

        _hw_real x_max_first_temp = u_curr_local + du_max[i]; //controled_system->acc_max[i];
		x_max_first[i] = (x_max_first[i] > u_max[i])? u_max[i]: x_max_first_temp;

		_hw_real x_min_first_temp = u_curr_local + (-du_max[i]);//controled_system->acc_min[i];        
		x_min_first[i] = (x_min_first[i] < u_min[i])? u_min[i]: x_min_first_temp;	
	}
}
// ---------------------------------------------------
void initializeBestLocalFitness(void){
#pragma HLS inline
	// Initialize best local fitness 
	for (int i = 0; i < _n_S; ++i) { 
		f_ind[i] = H_MAX;
	}
}
// ---------------------------------------------------
int detectGlobalMinimum(
	int iter
){
	int pos;

	// Global Minimum detection
    bestfitness[iter] = min_array(f_ind, &pos);
    
	// // Global minimum position
    // for (int i = 0; i < _Nu*_n_U; ++i) {
    //     global_min[i] = _y[pos][i];
    // }
    return pos;

}
// ---------------------------------------------------
void equalizeParticles(
	_hw_real _x[][_Nu*_n_U],
	_hw_real iteration
){
#pragma HLS inline
/*
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation
*/	
	_real x_tmp;
	if (iteration == 0){	
	// for (int i = 0; i < _n_S; ++i) {
	// 	for (int j = 0; j < _Nu; ++j) {
	// 		int idx = k*_n_U;
	// 		for (int k = 1; k < _n_U; ++k) {
	// 			_x[i][idx] = _x[i][j];
	// 			idx++;
    //         }
    //     }
    // }

	for (int j = 0; j < Nu; ++j) {
        for (int k = 1; k < Nc; ++k) {
            int idx = k*Nu+j;
            for (int i = 0; i < S; ++i) {
#pragma UNROLL
				x_tmp = _x[i][j];
                _x[i][idx] = x_tmp;
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
// #pragma HLS inline
	_hw_real local_value = value;
	_hw_real return_value;
	if (local_value > u_max[pos]) {
        return_value = u_max[pos];
    }else if (local_value < u_min[pos]){
		return_value = u_min[pos];
    }else{
		return_value = local_value;
	}
}

// ---------------------------------------------------
void detectInvalidParticles(
	int iter, 
	int best_pos,
	
	_hw_real _x[_n_S][_Nu*_n_U]
){
	for (int i = 0; i < _n_S; ++i) {
        if(_n_S != best_pos) {
            if((std::abs(fx[i]-bestfitness[iter])/bestfitness[iter]) < 0.0001){
//                int count = 0;  
                int count = _Nh*_n_U;
                for (int j = 0; j < _Nh*_n_U; ++j) {
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
	volatile _hw_real *uss,

	_hw_real _x[][_Nu*_n_U], 
	_hw_real _y[][_Nu*_n_U], 
	_hw_real _v[][_Nu*_n_U]
){
// INITIALIZATION OF PARTICLES WITH Delta u CONTRAINS FOR MPC
#pragma HLS inline
#pragma HLS ALLOCATION instances=rand_real limit=1 function
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation

	_hw_real x_ant[_n_U];
	memcpy_loop_rolled<_hw_real, _hw_real, _n_U>(x_ant, u_curr);

	for (int i = 0; i < _n_S; ++i) {
		for (int j = 0; j < _Nu; ++j) {
			int idx = j*_n_U;
			for (int k = 0; k < _n_U; ++k) {
#pragma HLS pipeline II=11 rewind
				_hw_real rand_tmp = rand_real();
		        _hw_real x_tmp = (x_ant[k] + (-du_max[k])) + ((_hw_real)2.0*du_max[k]) * rand_tmp; //random->read(); 
				_x[i][idx] = verifyControlConstrains(x_tmp, k);
                _y[i][idx] = uss[k];
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
    _hw_real _x[][_Nu*_n_U], 
    _hw_real _y[][_Nu*_n_U], 
    _hw_real _v[][_Nu*_n_U]
) {
// GENERAL PSO INITIALIZATION OF PARTICLES
	_hw_real d_max = 0.2;
    _hw_real d2_max = 0.4;
    int k = 0;
    for (int i = 0; i < _n_S; ++i) {
        //for (int j = 0; j < Nu/2; ++j) {
    	k = 0;
        for (int j = 0; j < _Nu; ++j) {
			_hw_real rand_tmp = rand_real();
            _hw_real x_tmp = (j == 0)? 
				u_min[k] + (u_max[k]-u_min[k]) * rand_tmp: //random->read();
				_x[i][j] = _x[i][j-1] - d_max + (d2_max) * rand_tmp; //random->read();
			_x[i][j] = x_tmp; 
			_y[i][j] = H_MAX;
            _v[i][j] = init_v;     
            k++; if (k>=_n_U) k=0;
        }
        valid_particle[i] = 1;
    }
}

// ---------------------------------------------------
void initializeLastBestKPSO(
	_hw_real volatile *last_best,
	
	_hw_real _x[][_Nu*_n_U]
){
#pragma HLS inline
	// Uses KPSO (best last position)
	for (int i = 0; i < _n_U; i++){
		_hw_real last_best_tmp = last_best[i];
		_hw_real x_local;
		if (last_best_tmp > x_max_first[i])
			x_local = x_max_first[i];
		else if (last_best_tmp < x_min_first[i])
			x_local = x_min_first[i];
		else
			x_local = last_best_tmp;
		_x[0][i] = x_local;
	}
	
	for (int i = _n_U; i < (_Nu*_n_U); i++) {
		_x[0][i] = last_best[i];			
	}
}
// ---------------------------------------------------
void initializeStableZero(
	_hw_real *uss, 
	_hw_real *u_curr, 
	int index,
	
	_hw_real _x[][_Nu*_n_U]
){
#pragma HLS inline
// Starts one particle with all Zeros for 
// stable response after equilibrium is reached
    int idx;
   	for (int k = 0; k < _n_U; ++k){
        for (int i = 0; i < _Nu; ++i) {
            idx = k*_Nu+i;
            _x[index][idx] = uss[0];
            if(i==0){
				_hw_real tmp_sub = _x[index][idx] - u_curr[k];
                if(tmp_sub > du_max[k]){
                    _x[index][idx] = du_max[k];
                }
                else if(tmp_sub < - du_max[k]){
                    _x[index][idx] = - du_max[k];
                }
            }
            else{
				_hw_real tmp_sub = _x[index][idx] - _x[index][idx];
                if(tmp_sub > du_max[k]){
                    _x[index][idx] = du_max[k];
                }
                else if(tmp_sub < - du_max[k]){
                    _x[index][idx] = - du_max[k];
                }
            }
            verifyControlConstrains(&_x[index][idx], k); 
        }
   	}
}
// ---------------------------------------------------
void evaluateFitnessAndDetectLocalBest(
	_hw_real _x[_n_S][_Nu*_n_U], 
	_hw_real _y[_n_S][_Nu*_n_U],

	_hw_real xref[_Nx*_Nu], 
	_hw_real uref[_n_U],
	_hw_real xss[_Nx]
){
#pragma HLS interface ap_bus  port=_x
#pragma HLS interface ap_bus  port=xref
#pragma HLS interface ap_fifo port=uref
#pragma HLS interface ap_fifo port=xss

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

	// Evaluates fitness and local detection
    loop_1: for(int i = 0; i < _n_S; i++) {
        // std::cout << std::endl;
#pragma HLS unroll factor=4 region skip_exit_check
        if(valid_particle[i] == 1){
            // std::cout << x[i][0] << " " << x[i][1] << " " << x[i][2] << " " << x[i][3];
            // fx[i] = rand_real(); 
			current_system.nmpc_cost_function(&_x[i][0], xref, uref, xss, uss);
            if (fx[i] < f_ind[i]) {
                loop_2: for (int j = 0; j < _Nu*_n_U; ++j) {
                	_y[i][j] = _x[i][j];
                }
                f_ind[i] = fx[i];           
            }
        }
    }
    // std::cout << std::endl;

}
// ---------------------------------------------------
void updateParticles(
	_hw_real global_min[],

    _hw_real _x[][_Nu*_n_U],
    _hw_real _y[][_Nu*_n_U],
    _hw_real _v[][_Nu*_n_U]
){
    _hw_real r1, r2;
    
    // Update Particles
    for (int i = 0; i < _n_S; ++i) {
        if(valid_particle[i] == 1){
            for (int j = 0; j < _Nu; ++j) {
                r1 = rand_real();
                r2 = rand_real();

                _v[i][j] = w*_v[i][j] + c1*r1*(_y[i][j]-_x[i][j]) + c2*r2*(global_min[j] - _x[i][j]);

                if (std::abs(_v[i][j]) > max_v) {
                    if (_v[i][j] > 0)
                        _v[i][j] = max_v;
                    else
                        _v[i][j] = -max_v;
                }
                _x[i][j] = _x[i][j] + _v[i][j];

                if (_x[i][j] > u_max[j]) {
                    _x[i][j] = u_max[j];
                }
                else if (_x[i][j] < u_min[j]){
                    _x[i][j] = u_min[j];
                }
            } // END FOR j
        } // End if
    } // END FOR i
}
// ---------------------------------------------------
void updateParticlesWithDuConstrains(
	_hw_real global_min[],

	_hw_real _x[_n_S][_Nu*_n_U], 
	_hw_real _y[_n_S][_Nu*_n_U], 
	_hw_real _v[_n_S][_Nu*_n_U]
){
#pragma HLS interface ap_bus  port=x
#pragma HLS interface ap_bus  port=y
#pragma HLS interface ap_bus  port=v

#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation

#pragma HLS ALLOCATION instances=rand_real limit=1 function

    _hw_real r1, r2;
    
    // Update Particles
    for (int i = 0; i < _n_S; ++i) {
        if(valid_particle[i] == 1){
    	 for (int j = 0; j < _Nu; ++j) {
    	 	for (int k = 0; k < _n_U; ++k){
    	 		int idx = k*_Nu+j;
                r1 = rand_real(); //random->read();
                r2 = rand_real(); //random->read();

	            _hw_real v_ant = _v[i][idx];
	            _v[i][idx] = w*_v[i][idx] + c1*r1*(_y[i][idx]-_x[i][idx]) + c2*r2*(global_min[idx] - _x[i][idx]);
	            if (std::abs(_v[i][idx]) > max_v) {
	                if (_v[i][idx] > 0)
	                    _v[i][idx] = max_v;
	                else
	                    _v[i][idx] = -max_v;
	            }
	            _x[i][idx] = _x[i][idx] + _v[i][idx];
                
	            if (j==0) {
	                if (_x[i][idx] > x_max_first[k]) {
	                    _x[i][idx] = x_max_first[k];
	                }
	                else if (_x[i][idx] < x_min_first[k]) {
	                    _x[i][idx] = x_min_first[k];
	                }
	            }
	            else {
	                if (std::abs(_x[i][idx]-_x[i][idx-1]) > du_max[k]) {
	                    if ((_x[i][idx]-_x[i][idx-1]) > 0)
	                        _x[i][idx] = _x[i][idx-1] + du_max[k];
	                    else
	                        _x[i][idx] = _x[i][idx-1] - du_max[k];
	                }
                }
                verifyControlConstrains(&_x[i][idx], k);
	            
            } // END FOR k
        } // END FOR j
       } // END IF
    } // END FOR i
}
};

#endif