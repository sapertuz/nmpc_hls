#include <float.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string> 
#include <sstream>

#include "pso.hpp"

// DEBUG
int itr;

PSO::PSO() {

	// Accessing config file

    // General PSO variables
    max_v = _max_v;		// Maximum velocity
	w0 = _w0; 				// initial weight	
	wf = _wf; 				// final weight
	w = w0;
	slope = (wf-w0)/maxiter;
	c1 = _c1; 				// cognitive coefficient
	c2 = _c2; 				// social coefficient
	threshold = _threshold;
	stop_criteria = _stop_criteria;
    
    // Controlled System Dependent Variables
    //this->controled_system = new ModelState();
    //std::cout << this->controled_system << std::endl;
    
    N = _N;
    Nu = _Nu;
    Nc = _n_U;
    for (_uchar i = 0; i < Nu; i++){
        x_min[i] = _u_min[i];
        x_max[i] = _u_max[i];
        du_max[i] = _du_max[i];
    }

    ini_v = max_v/10; // initial velocity

    number_of_active_particles = S;

    /* initialize random seed: */
    //srand (time(NULL));
#ifndef __arm__
	//sranddev();
#endif
}

PSO::~PSO() {
    std::cout<<"salir"<<std::endl;
}

int PSO::getS(){return S;}
int PSO::getKPSO(){return kpso;}
int PSO::getStableZero(){return stable_zero;}
int PSO::getMaxiter(){return maxiter;}
_real PSO::getMaxV(){return max_v;}
int PSO::getStopCriteria(){return stop_criteria;}
_real PSO::getC1(){return c1;}
_real PSO::getC2(){return c2;}

void PSO::initializeConstrains(
    _real u_curr[]
){
	// Initialize constrains based on current state
    
	for (int i = 0; i < Nc; ++i){
        x_max_first[i] = u_curr[i] + controled_system->acc_max[i];
        x_min_first[i] = u_curr[i] + controled_system->acc_min[i];
        
		if (x_max_first[i] > x_max[i])
		    x_max_first[i] = x_max[i];
		if (x_min_first[i] < x_min[i])
		    x_min_first[i] = x_min[i];	
	}
}

void PSO::verifyControlConstrains(
    _real value[], 
    int pos
){
    if (value[0] > x_max[pos]) {
        value[0] = x_max[pos];
    }
    else if (value[0] < x_min[pos]){
            value[0] = x_min[pos];
    }
}

// GENERAL PSO INITIALIZATION OF PARTICLES
void PSO::initializeParticles() {
    float d_max = 0.2;
    float d2_max = 0.4;
    for (int i = 0; i < S; ++i) {
        //for (int j = 0; j < Nu/2; ++j) {
        for (int j = 0; j < Nu; ++j) {
            if(j == 0){
                x[i][j] = x_min[j] + (x_max[j]-x_min[j]) * rand_real(); //random->read();
                y[i][j] = 1e20;
                v[i][j] = ini_v;
            }
            else{
                x[i][j] = x[i][j-1] - d_max + (d2_max) * rand_real(); //random->read();
                y[i][j] = 1e20;
                v[i][j] = ini_v;                
            }
            // x[i][j+2] = x[i][j];
            // y[i][j+2] = 1e20;
            // v[i][j+2] = ini_v;

        }
        valid_particle[i] = 1;
    }
}

// INITIALIZATION OF PARTICLES WITH Delta u CONTRAINS FOR MPC
void PSO::initializeParticlesWithDuConstrains(
    _real u_curr[], 
    _real uss[]) 
{
	for (int j = 0; j < Nu; ++j) {
		for (int k = 0; k < Nc; ++k) {
			int idx = k*Nu+j;
			for (int i = 0; i < S; ++i) {
		        if (j==0) {
		            //x[i][idx] = (u_curr[k] - du_max[k]) + (2*du_max[k]) * random->read(); //rand_real();
                    x[i][idx] = (u_curr[k] + controled_system->acc_min[k]) + (controled_system->acc_max[k]-controled_system->acc_min[k]) * rand_real(); //random->read(); 
		        }
		        else {
		            x[i][idx] = (x[i][idx-1] - du_max[k]) + (2*du_max[k]) * rand_real(); //random->read(); 
		        }
                verifyControlConstrains(&x[i][idx], k);
                // if (x[i][idx] > x_max[k]) {
                //     x[i][idx] = x_max[k];
                // }
                // else if (x[i][idx] < x_min[k]){
                //         x[i][idx] = x_min[k];
                // }

                y[i][idx] = uss[k];
		        v[i][idx] = ini_v;
                valid_particle[i] = 1;
		    }
		}
	}
}

void PSO::equalizeParticlesVant() {
    for (int j = 0; j < Nu; ++j) {
        for (int k = 1; k < Nc; ++k) {
            int idx = k*Nu+j;
            for (int i = 0; i < S; ++i) {
                x[i][idx] = x[i][j];
            }
        }
    }
}

void PSO::initializeLastBestKPSO(
    _real last_best[]
) {
	// Uses KPSO (best last position)
	for (int j = 0; j < Nc; ++j) {
	   	for (int i = 0; i < (Nu-1); ++i) {
	   		x[0][j*Nu+i] = last_best[j*Nu+i+1];
            if(i==0){
                if (x[0][j*Nu+i] > x_max_first[j])
                    x[0][j*Nu+i] = x_max_first[j];
                if (x[0][j*Nu+i] < x_min_first[j])
                    x[0][j*Nu+i] = x_min_first[j];
            }
                
	   	}
	   	x[0][(j+1)*Nu-1] = last_best[(j+1)*Nu-1];        
	}
}

void PSO::initializeStableZero(
    _real uss[], 
    _real u_curr[], 
    int index
){
    // Starts one particle with all Zeros for 
    // stable response after equilibrium is reached
    int idx;
   	for (int k = 0; k < Nc; ++k){
        for (int i = 0; i < Nu; ++i) {
            idx = k*Nu+i;
            x[index][idx] = uss[0];
            if(i==0){
                if((x[index][idx] - u_curr[k]) > du_max[k]){
                    x[index][idx] = du_max[k];
                }
                else if((x[index][idx] - u_curr[k]) < -du_max[k]){
                    x[index][idx] = -du_max[k];
                }
            }
            else{
                if((x[index][idx] - x[index][idx]) > du_max[k]){
                    x[index][idx] = du_max[k];
                }
                else if((x[index][idx] - x[index][idx]) < -du_max[k]){
                    x[index][idx] = -du_max[k];
                }
            }
            verifyControlConstrains(&x[index][idx], k); 
        }
   	}
}

void PSO::initializeBestLocalFitness(void){
	// Initialize best local fitness 
	for (int i = 0; i < S; ++i) { 
		f_ind[i] = FLT_MAX;
	}
}

void PSO::evaluateFitnessAndDetectLocalBest(
    _real x[][_Nc*_Nuu], 
    _real xref[], 
    _real uref[], 
    _real xss[], 
    _real uss[]
) {
	// Evaluates fitness and local detection
    for(int i = 0; i < S; i++) {
        // std::cout << std::endl;
        if(valid_particle[i] == 1){
            // std::cout << x[i][0] << " " << x[i][1] << " " << x[i][2] << " " << x[i][3];
            fx[i] = controled_system->nmpc_cost_function(&x[i][0], xref, uref, xss, uss);
            if (fx[i] < f_ind[i]) {
                for (int j = 0; j < Nu*Nc; ++j) {
                	y[i][j] = x[i][j];
                }
                f_ind[i] = fx[i];           
            }
        }
    }
    // std::cout << std::endl;
}

void PSO::detectInvalidParticles(
    int iter, 
    int best_pos
){
    for (int i = 0; i < S; ++i) {
        if(S != best_pos) {
            if((std::abs(fx[i]-bestfitness[iter])/bestfitness[iter]) < 0.0001){
//                int count = 0;  
                int count = N*Nc;
                for (int j = 0; j < N*Nc; ++j) {
                    if((std::abs(x[i][j]-global_min[j])/global_min[j]) < 0.0001){
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

int PSO::detectGlobalMinimum(
int iter
){
	int pos;

	// Global Minimum detection
    bestfitness[iter] = min_array(f_ind, &pos);
    
	// Global minimum position
    for (int i = 0; i < Nu*Nc; ++i) {
        global_min[i] = y[pos][i];
    }
    return pos;
}

void PSO::updateParticles(
    _real x[][_Nuu*_Nc], 
    _real y[][_Nuu*_Nc], 
    _real v[][_Nuu*_Nc]
){
    _real r1, r2;
    
    // Update Particles
    for (int i = 0; i < S; ++i) {
        if(valid_particle[i] == 1){
            for (int j = 0; j < Nu; ++j) {
                r1 = rand_real();
                r2 = rand_real();

                v[i][j] = w*v[i][j] + c1*r1*(y[i][j]-x[i][j]) + c2*r2*(global_min[j] - x[i][j]);

                if (std::abs(v[i][j]) > max_v) {
                    if (v[i][j] > 0)
                        v[i][j] = max_v;
                    else
                        v[i][j] = -max_v;
                }
                x[i][j] = x[i][j] + v[i][j];

                if (x[i][j] > x_max[j]) {
                    x[i][j] = x_max[j];
                }
                else if (x[i][j] < x_min[j]){
                    x[i][j] = x_min[j];
                }
            } // END FOR j
        } // End if
    } // END FOR i
}

void PSO::updateParticlesWithDuConstrains(
    _real x[][_Nuu*_Nc], 
    _real y[][_Nuu*_Nc], 
    _real v[][_Nuu*_Nc]
){
    _real r1, r2;
    
    // Update Particles
    for (int i = 0; i < S; ++i) {
        if(valid_particle[i] == 1){
    	 for (int j = 0; j < Nu; ++j) {
    	 	for (int k = 0; k < Nc; ++k){
    	 		int idx = k*Nu+j;
                r1 = rand_real(); //random->read();
                r2 = rand_real(); //random->read();

	            _real v_ant = v[i][idx];
	            v[i][idx] = w*v[i][idx] + c1*r1*(y[i][idx]-x[i][idx]) + c2*r2*(global_min[idx] - x[i][idx]);
	            if (std::abs(v[i][idx]) > max_v) {
	                if (v[i][idx] > 0)
	                    v[i][idx] = max_v;
	                else
	                    v[i][idx] = -max_v;
	            }
	            x[i][idx] = x[i][idx] + v[i][idx];
                
	            if (j==0) {
	                if (x[i][idx] > x_max_first[k]) {
	                    x[i][idx] = x_max_first[k];
	                }
	                else if (x[i][idx] < x_min_first[k]) {
	                    x[i][idx] = x_min_first[k];
	                }
	            }
	            else {
	                if (std::abs(x[i][idx]-x[i][idx-1]) > du_max[k]) {
	                    if ((x[i][idx]-x[i][idx-1]) > 0)
	                        x[i][idx] = x[i][idx-1] + du_max[k];
	                    else
	                        x[i][idx] = x[i][idx-1] - du_max[k];
	                }
                }
                verifyControlConstrains(&x[i][idx], k);
	            
            } // END FOR k
        } // END FOR j
       } // END IF
    } // END FOR i
}

int PSO::execute(   
    _real x_curr[], 
    _real u_curr[], 
    int iteration, 
    _real last_best[], 
    _real xref[], 
    _real uref[], 
    _real xss[], 
    _real uss[], 
    _real new_best[], 
    _real * J
){
    // Update sensor read
    this->controled_system->setState(x_curr);

    number_of_active_particles = S;
    itr = iteration;
    
	int best_pos;
    w0 = 0.9; // initial weight
	wf = 0.1; // final weight
	w = w0;
	slope = (wf-w0)/maxiter;
	ini_v = max_v/10;

	initializeConstrains(u_curr);

    initializeParticlesWithDuConstrains(u_curr, uss);

    if((iteration == 0) && (_n_U > 1)) {
        equalizeParticlesVant();
    }        

    if ((kpso == 1) && (iteration > 0)) {
        initializeLastBestKPSO(last_best);
    }
    if (stable_zero == 1) {
        initializeStableZero(uss, u_curr, 1);
        //initializeStableZero(x_max, u_curr, 2);
        //initializeStableZero(x_min, u_curr, 3);
    } 

    initializeBestLocalFitness();

	// ITERATIVE PROCESS
	int k = 0;  // index of iteration
	while (k < maxiter) {

        evaluateFitnessAndDetectLocalBest(x, xref, uref, xss, uss);     
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
    for (int i = 0; i < Nu*Nc; ++i) {
        new_best[i] = global_min[i];
    }

	*J = bestfitness[k-1];

	return k;
}

_real PSO::rand_real(){
	return (_real)rand()/(_real)RAND_MAX;
}

_real PSO::min_array(_real array[], int * pos) {
	//[bestfitness(k), p] = min(f_ind);
	_real min = array[0];
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
