#include <float.h>
#include <cmath>

#include "hls_pso.hpp"

// DEBUG
void hls_pso__init(HLS_PSO * myself) {

	// Accessing config file

    // General PSO variables
    myself->max_v = _max_v;		// Maximum velocity
	myself->w0 = _w0; 				// initial weight	
	myself->wf = _wf; 				// final weight
	myself->w = myself->w0;
	myself->slope = myself->slope_init;
	myself->c1 = _c1; 				// cognitive coefficient
	myself->c2 = _c2; 				// social coefficient
	myself->threshold = _threshold;
	myself->stop_criteria = _stop_criteria;
    
    // Controlled System Dependent Variables
    
    myself->N = _Nn;
    myself->Nu = _Nuu;
    myself->Nc = _Nc;
    for (int i = 0; i < _Nc; i++){
        myself->x_min[i] = _u_min[i];
        myself->x_max[i] = _u_max[i];
        myself->du_max[i] = _du_max[i];
    }

    myself->number_of_active_particles = _S;
}

void hls_pso__initializeConstrains(
    _hw_real u_curr[],
    _hw_real x_max_first[],
    _hw_real x_min_first[]
){
	// Initialize constrains based on current state
#pragma HLS inline

	for (int i = 0; i < _Nc; ++i){
        x_max_first[i] = u_curr[i] + _du_max[i]; //controled_system->acc_max[i];
        x_min_first[i] = u_curr[i] + (-_du_max[i]);//controled_system->acc_min[i];
        
		if (x_max_first[i] > _u_max[i])
		    x_max_first[i] = _u_max[i];
		if (x_min_first[i] < _u_min[i])
		    x_min_first[i] = _u_min[i];	
	}
}

void hls_pso__verifyControlConstrains(
    _hw_real value[],
    int pos
){

    if (value[0] > _u_max[pos]) {
        value[0] = _u_max[pos];
    }
    else if (value[0] < _u_min[pos]){
            value[0] = _u_min[pos];
    }
}

// GENERAL PSO INITIALIZATION OF PARTICLES
void hls_pso__initializeParticles(
    HLS_PSO * myself,

    _hw_real _x[][_Nu*_n_U], 
    _hw_real _y[][_Nu*_n_U], 
    _hw_real _v[][_Nu*_n_U]
) {
    float d_max = 0.2;
    float d2_max = 0.4;
    int k = 0;
    for (int i = 0; i < _S; ++i) {
        //for (int j = 0; j < Nu/2; ++j) {
    	k = 0;
        for (int j = 0; j < _Nu; ++j) {
            if(j == 0){
                _x[i][j] = _u_min[k] + (_u_max[k]-_u_min[k]) * hls_pso__rand_real(); //random->read();
                _y[i][j] = H_MAX;
                _v[i][j] = myself->ini_v;
            }
            else{
                _x[i][j] = _x[i][j-1] - d_max + (d2_max) * hls_pso__rand_real(); //random->read();
                _y[i][j] = H_MAX;
                _v[i][j] = myself->ini_v;                
            }
            k++; if (k>=_Nc) k=0;
        }
        myself->valid_particle[i] = 1;
    }
}

// INITIALIZATION OF PARTICLES WITH Delta u CONTRAINS FOR MPC
void hls_pso__initializeParticlesWithDuConstrains(
    HLS_PSO * myself,

    _hw_real u_curr[],
    _hw_real uss[],
    
    _hw_real _x[][_Nu*_n_U], 
    _hw_real _y[][_Nu*_n_U], 
    _hw_real _v[][_Nu*_n_U]
){
#pragma HLS inline

#pragma HLS ALLOCATION instances=rand_real limit=1 function

	for (int j = 0; j < _Nu; ++j) {
		for (int k = 0; k < _Nc; ++k) {
			int idx = k*_Nu+j;
			for (int i = 0; i < _S; ++i) {
		        if (j==0) {
                    //x[i][idx] = (u_curr[k] + controled_system->acc_min[k]) + (controled_system->acc_max[k]-controled_system->acc_min[k]) * rand_real(); //random->read();
                    _x[i][idx] = (u_curr[k] + (-_du_max[k])) + (_du_max[k]-(-_du_max[k])) * hls_pso__rand_real(); //random->read(); 
		        }
		        else {
		            _x[i][idx] = (_x[i][idx-1] - _du_max[k]) + (2*_du_max[k]) * hls_pso__rand_real(); //random->read(); 
		        }
                hls_pso__verifyControlConstrains(&_x[i][idx], k);

                _y[i][idx] = uss[k];
		        _v[i][idx] = myself->ini_v;
                myself->valid_particle[i] = 1;
		    }
		}
	}
}

void hls_pso__equalizeParticlesVant(
    _hw_real _x[][_Nu*_n_U]
) {
#pragma HLS inline
    for (int j = 0; j < _Nu; ++j) {
        for (int k = 1; k < _Nc; ++k) {
            int idx = k*_Nu+j;
            for (int i = 0; i < _S; ++i) {
                _x[i][idx] = _x[i][j];
            }
        }
    }
}

void hls_pso__initializeLastBestKPSO(
    HLS_PSO * myself,

    _hw_real last_best[],
    
    _hw_real _x[][_Nu*_n_U]
) {
#pragma HLS interface ap_bus  port=last_best

#pragma HLS inline
	// Uses KPSO (best last position)
	for (int j = 0; j < _Nc; ++j) {
	   	for (int i = 0; i < (_Nu-1); ++i) {
	   		_x[0][j*_Nu+i] = last_best[j*_Nu+i+1];
            if(i==0){
                if (_x[0][j*_Nu+i] > myself->x_max_first[j])
                    _x[0][j*_Nu+i] = myself->x_max_first[j];
                if (_x[0][j*_Nu+i] < myself->x_min_first[j])
                    _x[0][j*_Nu+i] = myself->x_min_first[j];
            }
	   	}
	   	_x[0][(j+1)*_Nu-1] = last_best[(j+1)*_Nu-1];        
	}
}

void hls_pso__initializeStableZero(
    _hw_real uss[],
    _hw_real u_curr[],
    int index,
    
    _hw_real _x[][_Nu*_n_U]
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
                if((_x[index][idx] - u_curr[k]) > _du_max[k]){
                    _x[index][idx] = _du_max[k];
                }
                else if((_x[index][idx] - u_curr[k]) < -_du_max[k]){
                    _x[index][idx] = -_du_max[k];
                }
            }
            else{
                if((_x[index][idx] - _x[index][idx]) > _du_max[k]){
                    _x[index][idx] = _du_max[k];
                }
                else if((_x[index][idx] - _x[index][idx]) < -_du_max[k]){
                    _x[index][idx] = -_du_max[k];
                }
            }
            hls_pso__verifyControlConstrains(&_x[index][idx], k); 
        }
   	}
}

void hls_pso__initializeBestLocalFitness(
    _hw_real f_ind[]
){
#pragma HLS inline
	// Initialize best local fitness 
	for (int i = 0; i < _S; ++i) { 
		f_ind[i] = H_MAX;
	}
}

void hls_pso__evaluateFitnessAndDetectLocalBest(
    HLS_PSO * myself,
    
    _hw_real _x[][_Nc*_Nuu],
    _hw_real _y[_S][_Nu*_n_U], 

    _hw_real xref[],
    _hw_real uref[],
    _hw_real xss[],
    _hw_real uss[]
) {
#pragma HLS interface ap_bus  port=_x
#pragma HLS interface ap_bus  port=xref
#pragma HLS interface ap_fifo port=uref
#pragma HLS interface ap_fifo port=xss
#pragma HLS interface ap_fifo port=uss

	// Evaluates fitness and local detection
    loop_1: for(int i = 0; i < _S; i++) {
        // std::cout << std::endl;
#pragma HLS unroll factor=4 region skip_exit_check
        if(myself->valid_particle[i] == 1){
            // std::cout << x[i][0] << " " << x[i][1] << " " << x[i][2] << " " << x[i][3];
            myself->fx[i] = hls_pso__rand_real(); // controled_system->nmpc_cost_function(&_x[i][0], xref, uref, xss, uss);
            if (myself->fx[i] < myself->f_ind[i]) {
                loop_2: for (int j = 0; j < _Nu*_Nc; ++j) {
                	_y[i][j] = _x[i][j];
                }
                myself->f_ind[i] = myself->fx[i];           
            }
        }
    }
    // std::cout << std::endl;
}

void hls_pso__detectInvalidParticles(
    HLS_PSO * myself,

    int iter, 
    int best_pos,

    _hw_real bestfitness[_maxiter],
    _hw_real global_min[_Nu*_n_U],
    _hw_real _x[_S][_Nu*_n_U]
){
    for (int i = 0; i < _S; ++i) {
        if(_S != best_pos) {
            if((std::abs(myself->fx[i]-bestfitness[iter])/bestfitness[iter]) < 0.0001){
//                int count = 0;  
                int count = _N*_Nc;
                for (int j = 0; j < _N*_Nc; ++j) {
                    if((std::abs(_x[i][j]-global_min[j])/global_min[j]) < 0.0001){
                        count--;
                    }
                }
                if((myself->valid_particle[i] == 1) && (count == 0)){
                    myself->valid_particle[i] = 0;
                    myself->number_of_active_particles--;
                }
            }
        }
    }
}

int hls_pso__detectGlobalMinimum(
	HLS_PSO * myself,
    int iter,
    
	_hw_real _y[][_Nu*_n_U],
    _hw_real bestfitness[],
    _hw_real global_min[]
){
	int pos;

	// Global Minimum detection
    bestfitness[iter] = hls_pso__min_array(myself->f_ind, &pos);
    
	// Global minimum position
    for (int i = 0; i < _Nu*_Nc; ++i) {
        global_min[i] = _y[pos][i];
    }
    return pos;
}

void hls_pso__updateParticles(
	HLS_PSO * myself,

    _hw_real global_min[],

    _hw_real _x[][_Nuu*_Nc],
    _hw_real _y[][_Nuu*_Nc],
    _hw_real _v[][_Nuu*_Nc]
){
    _hw_real r1, r2;
    
    // Update Particles
    for (int i = 0; i < _S; ++i) {
        if(myself->valid_particle[i] == 1){
            for (int j = 0; j < _Nu; ++j) {
                r1 = hls_pso__rand_real();
                r2 = hls_pso__rand_real();

                _v[i][j] = myself->w*_v[i][j] + myself->c1*r1*(_y[i][j]-_x[i][j]) + myself->c2*r2*(global_min[j] - _x[i][j]);

                if (std::abs(_v[i][j]) > myself->max_v) {
                    if (_v[i][j] > 0)
                        _v[i][j] = myself->max_v;
                    else
                        _v[i][j] = -myself->max_v;
                }
                _x[i][j] = _x[i][j] + _v[i][j];

                if (_x[i][j] > myself->x_max[j]) {
                    _x[i][j] = myself->x_max[j];
                }
                else if (_x[i][j] < myself->x_min[j]){
                    _x[i][j] = myself->x_min[j];
                }
            } // END FOR j
        } // End if
    } // END FOR i
}

void hls_pso__updateParticlesWithDuConstrains(
    HLS_PSO * myself,

    _hw_real global_min[],
    
    _hw_real _x[][_Nuu*_Nc],
    _hw_real _y[][_Nuu*_Nc],
    _hw_real _v[][_Nuu*_Nc]
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
        if(myself->valid_particle[i] == 1){
    	 for (int j = 0; j < _Nu; ++j) {
    	 	for (int k = 0; k < _Nc; ++k){
    	 		int idx = k*_Nu+j;
                r1 = hls_pso__rand_real(); //random->read();
                r2 = hls_pso__rand_real(); //random->read();

	            _hw_real v_ant = _v[i][idx];
	            _v[i][idx] = myself->w*_v[i][idx] + myself->c1*r1*(_y[i][idx]-_x[i][idx]) + myself->c2*r2*(global_min[idx] - _x[i][idx]);
	            if (std::abs(_v[i][idx]) > myself->max_v) {
	                if (_v[i][idx] > 0)
	                    _v[i][idx] = myself->max_v;
	                else
	                    _v[i][idx] = -myself->max_v;
	            }
	            _x[i][idx] = _x[i][idx] + _v[i][idx];
                
	            if (j==0) {
	                if (_x[i][idx] > myself->x_max_first[k]) {
	                    _x[i][idx] = myself->x_max_first[k];
	                }
	                else if (_x[i][idx] < myself->x_min_first[k]) {
	                    _x[i][idx] = myself->x_min_first[k];
	                }
	            }
	            else {
	                if (std::abs(_x[i][idx]-_x[i][idx-1]) > myself->du_max[k]) {
	                    if ((_x[i][idx]-_x[i][idx-1]) > 0)
	                        _x[i][idx] = _x[i][idx-1] + myself->du_max[k];
	                    else
	                        _x[i][idx] = _x[i][idx-1] - myself->du_max[k];
	                }
                }
                hls_pso__verifyControlConstrains(&_x[i][idx], k);
	            
            } // END FOR k
        } // END FOR j
       } // END IF
    } // END FOR i
}

int hls_pso__execute(
    HLS_PSO * myself,
    _hw_real x_curr[],
    _hw_real u_curr[],
    int iteration, 
    _hw_real last_best[],
    _hw_real xref[],
    _hw_real uref[],
    _hw_real xss[],
    _hw_real uss[],
    _hw_real new_best[],
    _hw_real * J
){
#pragma HLS ALLOCATION instances=hmul limit=1 operation
#pragma HLS ALLOCATION instances=hadd limit=1 operation
#pragma HLS ALLOCATION instances=hsub limit=1 operation
    // Update sensor read
    //static System * controled_system = new ModelState();
    //controled_system->setState(x_curr);

	// Particles Variables
	_hw_real x[_S][_Nu*_n_U];
	_hw_real y[_S][_Nu*_n_U];
	_hw_real v[_S][_Nu*_n_U];
	_hw_real bestfitness[_maxiter];
	_hw_real global_min[_Nu*_n_U];
    myself->number_of_active_particles = _S;
    
	int best_pos;

	hls_pso__initializeConstrains(u_curr);
    hls_pso__initializeParticlesWithDuConstrains(u_curr, uss, x, y, v);

    if((iteration == 0) && (_n_U > 1)) {
        hls_pso__equalizeParticles(x);
    }        

//    if ((kpso == 1) && (iteration > 0)) {
    hls_pso__initializeLastBestKPSO(last_best, x);
//    }
    if (myself->stable_zero == 1) {
        hls_pso__initializeStableZero(uss, u_curr, 1, x);
    } 

    hls_pso__initializeBestLocalFitness();
	// ITERATIVE PROCESS
	int k = 0;  // index of iteration

	while (k < myself->maxiter) {
        hls_pso__evaluateFitnessAndDetectLocalBest(x, y, xref, uref, xss, uss);     
        best_pos = hls_pso__detectGlobalMinimum(k);
        hls_pso__updateParticlesWithDuConstrains(x, y, v);    

        //detectInvalidParticles(k, best_pos);

		// STOP CRITERIA
		if(myself->stop_criteria) {
			if (k > myself->stop_criteria) {
			    if ((bestfitness[k-myself->stop_criteria]-bestfitness[k]) < myself->threshold)
			        break;
			}
		}
        if(myself->number_of_active_particles < 2){
            break;
        }
		k = k + 1;
		myself->w = myself->w + myself->slope;
	}
    
	// Return best value
    for (int i = 0; i < _Nu*_Nc; ++i) {
        new_best[i] = global_min[i];
    }
    
    *J = bestfitness[myself->maxiter-1];

	return myself->maxiter;
}

_hw_real hls_pso__rand_real(){
#ifndef __SYNTHESIS__
	return (_hw_real)rand()/(_hw_real)RAND_MAX;
#else
#pragma HLS ALLOCATION instances=sub limit=1 operation
#pragma HLS ALLOCATION instances=add limit=1 operation
#pragma HLS ALLOCATION instances=urem limit=1 operation
#pragma HLS ALLOCATION instances=udiv limit=1 operation
#pragma HLS ALLOCATION instances=srem limit=1 operation
#pragma HLS ALLOCATION instances=sdiv limit=1 operation
#pragma HLS ALLOCATION instances=mul limit=1 operation
    rnd_seed = (a * rnd_seed + c) % m;
    int tmp = rnd_seed%PSEUDORAND_MAX - PSEUDORAND_MAX_H;
    _hw_real return_value = (_hw_real)tmp*(_hw_real)PSEUDORAND_MAX_UNIT;
    return return_value;
#endif
}

_hw_real hls_pso__min_array(_hw_real array[], int * pos) {
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
