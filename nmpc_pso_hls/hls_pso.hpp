#ifndef HLS_PSO_HPP
#define HLS_PSO_HPP

//#define DEBUG_HLS

//#include "config.hpp"
#include "hls_system.hpp"

// #if _Parametrization > 0
//     #define _Nn     _Parametrization
//     #define _Nuu    _N
//     #define _Nc     1
// #else
//     #define _Nn     _N
//     #define _Nuu    _Nu
//     #define _Nc     _n_U
// #endif 

template<
    class _hw_real,
	class _system_t,
	class _randCore_t,
    unsigned _N,
    unsigned _Nx,
    unsigned _n_U,
    unsigned _Nu
>class PSO{
protected:
	// Pseudo random generator variables
	const int a = (1103515245);
	const int c = (12345);
	const int m = (1<<31);
	short rnd_seed = 9727;

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

	const _hw_real x_min[_Nn];
	const _hw_real x_max[_Nn];
	const _hw_real du_max[_n_U];
	
	_hw_real x_max_first[_n_U];
	_hw_real x_min_first[_n_U];

	const _hw_real init_v; // initial velocity

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
	
	// System 
	_system_t current_system;
	_randCore_t randGen;

public:
// Init Execute
constexpr PSO(
	// PSO Configuration
	const int 	 	_S,
	const int 	 	_stable_zero,
	const int 	 	_maxiter,
	const _hw_real  _max_v,
	const _hw_real  _w0,
	const _hw_real 	_wf,
	const _hw_real 	_slope,
	const _hw_real 	_c1,
	const _hw_real 	_c2,
	const _hw_real  _slope_init,
	const _hw_real	_x_min[_Nn],
	const _hw_real	_x_max[_Nn],
	const _hw_real	_du_max[_n_U],
	_system_t _current_system,
	_randCore_t _randGen
	) : S(_S), maxiter(_maxiter), max_v(_max_v), w0(_w0), wf(_wf), 
		slope(_slope), c1(_c1), c2(_c2), slope_init((_wf-_w0)/_maxiter),
		stable_zero(_stable_zero), init_v(_max_v*0.01),
		x_min(_x_min), x_max(_x_max), du_max(_du_max),
		current_system(_current_system), randGen(_randGen)
{
}

// ---------------------------------------------------
int execute(
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
){
    // Update sensor read
    //static System * controled_system = new ModelState();
    //controled_system->setState(x_curr);

	// Particles Variables
	_hw_real x[_S][_Nu*_n_U];
	_hw_real y[_S][_Nu*_n_U];
	_hw_real v[_S][_Nu*_n_U];
	_hw_real bestfitness[_maxiter];
	_hw_real global_min[_Nu*_n_U];
    number_of_active_particles = _S;
    
	int best_pos;

	initializeConstrains(u_curr);
    initializeParticlesWithDuConstrains(u_curr, uss, x, y, v);

    if((iteration == 0) && (_n_U > 1)) {
        equalizeParticles(x);
    }        

//    if ((kpso == 1) && (iteration > 0)) {
    initializeLastBestKPSO(last_best, x);
//    }
    if (stable_zero == 1) {
        initializeStableZero(uss, u_curr, 1, x);
    } 

    initializeBestLocalFitness();
	// ITERATIVE PROCESS
	int k = 0;  // index of iteration

	while (k < maxiter) {
        evaluateFitnessAndDetectLocalBest(x, y, xref, uref, xss, uss);     
        best_pos = detectGlobalMinimum(k);
        updateParticlesWithDuConstrains(x, y, v);    

        //detectInvalidParticles(k, best_pos);

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
		k = k + 1;
		w = w + slope;
	}
    
	// Return best value
    for (int i = 0; i < _Nu*_Nc; ++i) {
        new_best[i] = global_min[i];
    }
    
    *J = bestfitness[maxiter-1];

	return maxiter;
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
	_hw_real array[_S], 
	int * pos
){
	//[bestfitness(k), p] = min(f_ind);
#pragma HLS inline
	_hw_real min = array[0];
	int pos_temp = 0;
	for (int i = 1; i < _S; ++i) {
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
	_hw_real u_curr[_n_U]
){
	// Initialize constrains based on current state
#pragma HLS inline
	for (int i = 0; i < _Nc; ++i){
        x_max_first[i] = u_curr[i] + du_max[i]; //controled_system->acc_max[i];
        x_min_first[i] = u_curr[i] + (-du_max[i]);//controled_system->acc_min[i];        
		if (x_max_first[i] > u_max[i])
		    x_max_first[i] = u_max[i];
		if (x_min_first[i] < u_min[i])
		    x_min_first[i] = u_min[i];	
	}
	number_of_active_particles = S;
}
// ---------------------------------------------------
void initializeBestLocalFitness(void){
#pragma HLS inline
	// Initialize best local fitness 
	for (int i = 0; i < _S; ++i) { 
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
    
	// Global minimum position
    for (int i = 0; i < _Nu*_Nc; ++i) {
        global_min[i] = _y[pos][i];
    }
    return pos;

}
// ---------------------------------------------------
void equalizeParticles(
	_hw_real _x[_S][_Nu*_n_U]
){
#pragma HLS inline
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation
    for (int j = 0; j < _Nu; ++j) {
        for (int k = 1; k < _Nc; ++k) {
            int idx = k*_Nu+j;
            for (int i = 0; i < _S; ++i) {
                _x[i][idx] = _x[i][j];
            }
        }
    }

}
// ---------------------------------------------------
void verifyControlConstrains(
	_hw_real value[1], 
	int pos
){
#pragma HLS inline
	if (value[0] > u_max[pos]) {
        value[0] = u_max[pos];
    }
    else if (value[0] < u_min[pos]){
            value[0] = u_min[pos];
    }

}
// ---------------------------------------------------
void detectInvalidParticles(
	int iter, 
	int best_pos,
	
	_hw_real _x[_S][_Nu*_n_U]
){
	for (int i = 0; i < _S; ++i) {
        if(_S != best_pos) {
            if((std::abs(fx[i]-bestfitness[iter])/bestfitness[iter]) < 0.0001){
//                int count = 0;  
                int count = _N*_Nc;
                for (int j = 0; j < _N*_Nc; ++j) {
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
	_hw_real u_curr[_n_U], 
	_hw_real uss[_n_U],

	_hw_real _x[_S][_Nu*_n_U], 
	_hw_real _y[_S][_Nu*_n_U], 
	_hw_real _v[_S][_Nu*_n_U]
){
// INITIALIZATION OF PARTICLES WITH Delta u CONTRAINS FOR MPC
#pragma HLS inline
#pragma HLS ALLOCATION instances=rand_real limit=1 function
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation

	for (int j = 0; j < _Nu; ++j) {
		for (int k = 0; k < _Nc; ++k) {
			int idx = k*_Nu+j;
			for (int i = 0; i < _S; ++i) {
		        if (j==0) {
                    //x[i][idx] = (u_curr[k] + controled_system->acc_min[k]) + (controled_system->acc_max[k]-controled_system->acc_min[k]) * rand_real(); //random->read();
                    _x[i][idx] = (u_curr[k] + (-du_max[k])) + (du_max[k]-(-du_max[k])) * rand_real(); //random->read(); 
		        }
		        else {
		            _x[i][idx] = (_x[i][idx-1] - du_max[k]) + (2*du_max[k]) * rand_real(); //random->read(); 
		        }
                verifyControlConstrains(&_x[i][idx], k);

                _y[i][idx] = uss[k];
		        _v[i][idx] = init_v;
                valid_particle[i] = 1;
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
	float d_max = 0.2;
    float d2_max = 0.4;
    int k = 0;
    for (int i = 0; i < _S; ++i) {
        //for (int j = 0; j < Nu/2; ++j) {
    	k = 0;
        for (int j = 0; j < _Nu; ++j) {
            if(j == 0){
                _x[i][j] = _u_min[k] + (u_max[k]-u_min[k]) * rand_real(); //random->read();
                _y[i][j] = H_MAX;
                _v[i][j] = init_v;
            }
            else{
                _x[i][j] = _x[i][j-1] - d_max + (d2_max) * rand_real(); //random->read();
                _y[i][j] = H_MAX;
                _v[i][j] = init_v;                
            }
            k++; if (k>=_Nc) k=0;
        }
        valid_particle[i] = 1;
    }
}

// ---------------------------------------------------
void initializeLastBestKPSO(
	_hw_real last_best[_Nu*_n_U],
	
	_hw_real _x[_S][_Nu*_n_U]
){
#pragma HLS inline
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation
	// Uses KPSO (best last position)
	for (int j = 0; j < _Nc; ++j) {
	   	for (int i = 0; i < (_Nu-1); ++i) {
	   		_x[0][j*_Nu+i] = last_best[j*_Nu+i+1];
            if(i==0){
                if (_x[0][j*_Nu+i] > x_max_first[j])
                    _x[0][j*_Nu+i] = x_max_first[j];
                if (_x[0][j*_Nu+i] < x_min_first[j])
                    _x[0][j*_Nu+i] = x_min_first[j];
            }
	   	}
	   	_x[0][(j+1)*_Nu-1] = last_best[(j+1)*_Nu-1];        
	}
}
// ---------------------------------------------------
void initializeStableZero(
	_hw_real uss[_n_U], 
	_hw_real u_curr[_n_U], 
	int index,
	
	_hw_real _x[_S][_Nu*_n_U]
){
#pragma HLS inline
// Starts one particle with all Zeros for 
// stable response after equilibrium is reached
    int idx;
   	for (int k = 0; k < _Nc; ++k){
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
	_hw_real _x[_S][_Nu*_n_U], 
	_hw_real _y[_S][_Nu*_n_U],

	_hw_real xref[_Nx*_Nu], 
	_hw_real uref[_n_U],
	_hw_real xss[_Nx], 
	_hw_real uss[_n_U]
){
#pragma HLS interface ap_bus  port=_x
#pragma HLS interface ap_bus  port=xref
#pragma HLS interface ap_fifo port=uref
#pragma HLS interface ap_fifo port=xss
#pragma HLS interface ap_fifo port=uss
	typedef System<_real, _N, _Nx, _n_U, _Nu> T_system;
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
    loop_1: for(int i = 0; i < _S; i++) {
        // std::cout << std::endl;
#pragma HLS unroll factor=4 region skip_exit_check
        if(valid_particle[i] == 1){
            // std::cout << x[i][0] << " " << x[i][1] << " " << x[i][2] << " " << x[i][3];
            // fx[i] = rand_real(); 
			controled_system->nmpc_cost_function(&_x[i][0], xref, uref, xss, uss);
            if (fx[i] < f_ind[i]) {
                loop_2: for (int j = 0; j < _Nu*_Nc; ++j) {
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

    _hw_real _x[][_Nuu*_Nc],
    _hw_real _y[][_Nuu*_Nc],
    _hw_real _v[][_Nuu*_Nc]
){
    _hw_real r1, r2;
    
    // Update Particles
    for (int i = 0; i < _S; ++i) {
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

                if (_x[i][j] > x_max[j]) {
                    _x[i][j] = x_max[j];
                }
                else if (_x[i][j] < x_min[j]){
                    _x[i][j] = x_min[j];
                }
            } // END FOR j
        } // End if
    } // END FOR i
}
// ---------------------------------------------------
void updateParticlesWithDuConstrains(
	_hw_real global_min[],

	_hw_real _x[_S][_Nu*_n_U], 
	_hw_real _y[_S][_Nu*_n_U], 
	_hw_real _v[_S][_Nu*_n_U]
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
    for (int i = 0; i < _S; ++i) {
        if(valid_particle[i] == 1){
    	 for (int j = 0; j < _Nu; ++j) {
    	 	for (int k = 0; k < _Nc; ++k){
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