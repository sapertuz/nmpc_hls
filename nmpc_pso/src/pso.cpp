 #include "pso.hpp"
#include <float.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "aux_functions.hpp"
#include <string> 
#include <sstream>
#include "read_from_file.hpp"

// DEBUG
int itr;

PSO::PSO(System * controled_system, std::string config_file) {

	// Accessing config file
	std::ifstream sim_config;
    sim_config.open(config_file, std::ios::in);
    
    if(!sim_config.is_open())
    	throw std::runtime_error("PSO config file not found.");

    // General PSO variables
    S = read_int(&sim_config, "S");				// Prediction Horizon

	kpso = read_int(&sim_config, "KPSO");
	stable_zero = read_int(&sim_config, "stable_zero");;

	maxiter = read_int(&sim_config, "maxiter"); 	// Maximum number of iterations
	max_v = read_real(&sim_config, "max_v");		// Maximum velocity
	w0 = read_real(&sim_config, "w0"); 				// initial weight	
	wf = read_real(&sim_config, "wf"); 				// final weight
	w = w0;
	slope = (wf-w0)/maxiter;
	c1 = read_real(&sim_config, "c1"); 				// cognitive coefficient
	c2 = read_real(&sim_config, "c2"); 				// social coefficient
	threshold = read_real(&sim_config, "threshold");
	stop_criteria = read_int(&sim_config, "stop_criteria");
    
    // Controlled System Dependent Variables
    this->controled_system = controled_system;
    if(controled_system->getParametrized()) {
        N = controled_system->getParametrized();
        Nu = N;
        Nc = 1;

        x_min = controled_system->get_pmin();
        x_max = controled_system->get_pmax();
        max_v = 0.1;
        
        print_array("pmin", x_min, 4, 0);
        print_array("pmax", x_max, 4, 0);
        
        //du_max =  controled_system->getDU_max();
        previous_control = (_real *) malloc(controled_system->getParametrized()*sizeof(_real));
    }
    else{
        N = controled_system->getN();
        Nu = controled_system->getNu();
        Nc = controled_system->getn_U();

        x_min = controled_system->getU_min();
        x_max = controled_system->getU_max();
        du_max =  controled_system->getDU_max();        
    }

    if((x = alloc_matrix(S, Nu*Nc)) == NULL) {throw(1);}
    if((y = alloc_matrix(S, Nu*Nc)) == NULL) {throw(1);}
    if((v = alloc_matrix(S, Nu*Nc)) == NULL) {throw(1);}

    ini_v = max_v/10; // initial velocity

	x_max_first = (_real *) malloc(Nc*sizeof(_real));   
	x_min_first = (_real *) malloc(Nc*sizeof(_real));   

	f_ind = (_real *) malloc(S*sizeof(_real));
	fx    = (_real *) malloc(S*sizeof(_real));
	bestfitness = (_real *) malloc(maxiter*sizeof(_real));
	global_min = (_real *) malloc(Nu*Nc*sizeof(_real));

	sim_config.close();

    u_from_parameters = (_real *) malloc(controled_system->getN()*controled_system->getn_U()*sizeof(_real));

    valid_particle = (int *) malloc(S*sizeof(int));
    number_of_active_particles = S;

    /* initialize random seed: */
    //srand (time(NULL));
#ifndef __arm__
	//sranddev();
#endif
}

PSO::~PSO() {
	if(x != NULL) free_matrix(x, S);
	if(y != NULL) free_matrix(y, S);
	if(v != NULL) free_matrix(v, S);
	if(f_ind != NULL) free(f_ind);
	if(fx != NULL) free(fx);
	if(bestfitness != NULL) free(bestfitness);
	if(global_min != NULL) free(global_min);
    if(valid_particle != NULL) free(valid_particle);
	//delete(exporter);
}

int PSO::getS(){return S;}
int PSO::getKPSO(){return kpso;}
int PSO::getStableZero(){return stable_zero;}
int PSO::getMaxiter(){return maxiter;}
_real PSO::getMaxV(){return max_v;}
int PSO::getStopCriteria(){return stop_criteria;}
_real PSO::getC1(){return c1;}
_real PSO::getC2(){return c2;}

_real * PSO::get_x_max_first(){return x_max_first;}
_real * PSO::get_x_min_first(){return x_min_first;}
_real * PSO::get_f_ind(){return f_ind;}
_real * PSO::get_global_min(){return global_min;}

void PSO::initializeConstrains(_real * u_curr){
	// Initialize constrains based on current state
    
	for (int i = 0; i < Nc; ++i)
	{
		//x_max_first[i] = u_curr[i] + du_max[i];
		//x_min_first[i] = u_curr[i] - du_max[i];
        
        x_max_first[i] = u_curr[i] + controled_system->acc_max[i];
        x_min_first[i] = u_curr[i] + controled_system->acc_min[i];

        //std::cout << "x_max_first: " << x_max_first[i] << "  x_min_first: " << x_min_first[i] << std::endl;

		if (x_max_first[i] > x_max[i])
		    x_max_first[i] = x_max[i];
		if (x_min_first[i] < x_min[i])
		    x_min_first[i] = x_min[i];	
	}
}

void PSO::verifyControlConstrains(_real * value, int pos){
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
void PSO::initializeParticlesWithDuConstrains(_real * u_curr, _real * uss) {
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

void PSO::initializeLastBestKPSOParameters(_real * last_best) {
    // Uses KPSO (best last position)
    for (int j = 0; j < N; ++j) {
            x[0][j] = last_best[j]; 
            y[0][j] = last_best[j]; 
            v[0][j] = 0.0;
    }   
}

void PSO::initializeLastBestKPSO(_real * last_best) {
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

void PSO::initializeStableZeroParameters(_real * uss, int index){
    // Starts one particle the initial conditions passed
    for (int i = 0; i < Nu; ++i) {
        x[index][i] = uss[i];
        y[index][i] = uss[i];
        v[index][i] = 0.0;

    }
}

void PSO::initializeStableZero(_real * uss, _real * u_curr, int index){
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
    //print_array("fx_init", f_ind, S, 0);
}

void PSO::evaluateFitnessAndDetectLocalBest(_real ** x, _real * xref, _real * uref, _real * xss, _real * uss) {
	// Evaluates fitness and local detection
    for(int i = 0; i < S; i++) {
        if(valid_particle[i] == 1){
            fx[i] = controled_system->nmpc_cost_function(&x[i][0], xref, uref, xss, uss);
            if (fx[i] < f_ind[i]) {
            	// printf("fx[%d] < f_ind[%d] | %f < %f \n", i, i, fx[i], f_ind[i]);
                for (int j = 0; j < Nu*Nc; ++j) {
                	y[i][j] = x[i][j];
                }
                f_ind[i] = fx[i];           
            }
        }
    }
//    if((itr > 42) && (itr < 46)) {
//        print_array("fx   ", fx, S, 1);
//        //print_array("f_ind", f_ind, S, 1);
//    }
}

void PSO::evaluateFitnessAndDetectLocalBestParameters(_real ** x, _real * xref, _real * uref, _real * xss, _real * uss) {
    // Evaluates fitness and local detection
    for(int i = 0; i < S; i++) {
        if(valid_particle[i] == 1){
            
            controled_system->control_from_parameters(&x[i][0], previous_control, u_from_parameters, controled_system->getN());
            fx[i] = controled_system->nmpc_cost_function(u_from_parameters, xref, uref, xss, uss);
#ifdef DEBUG
            if(isinf(fx[i]) || isnan(fx[i])) {
                printf("State: ");
                for (int j=0; j < controled_system->getNx(); j++) {
                    printf("%f ", controled_system->getState(j));
                }
                printf("\n");
                printf("fx[%d]: %f\n", i, fx[i]);
                print_array("u_prev", previous_control, 4, 0);
                print_array("x", &x[i][0], 4, 0);
                print_array("u1", &u_from_parameters[0], 50, 0);
                print_array("u2", &u_from_parameters[50], 50, 0);
                print_array("u3", &u_from_parameters[100], 50, 0);
                print_array("u4", &u_from_parameters[150], 50, 0);
            }
#endif
            if (fx[i] < f_ind[i]) {
//                printf("fx[%d] < f_ind[%d] | %f < %f \n", i, i, fx[i], f_ind[i]);
                for (int j = 0; j < Nu*Nc; ++j) {
                    y[i][j] = x[i][j];
                }
                f_ind[i] = fx[i];           
            }
        }
    }
    //print_array("fx", fx, S, 1);
}

void PSO::detectInvalidParticles(int iter, int best_pos){
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
#ifdef DEBUG
    print_array_int("Valid", valid_particle, S);
#endif
}

int PSO::detectGlobalMinimum(int iter) {
	int pos;

	// Global Minimum detection
    bestfitness[iter] = min_array(f_ind, &pos, S);
    
    //std::cout << "Best POS: " << pos << " value: " << f_ind[pos] << std::endl;
#ifdef DEBUG
    for (int i = 0; i < S; i++) {
        if(valid_particle[i] != 0)
            valid_particle[i] = 1;
    }
    valid_particle[pos] = 2;
#endif
    
	// Global minimum position
    for (int i = 0; i < Nu*Nc; ++i) {
        global_min[i] = y[pos][i];
    }
    //print_array("global_min", global_min, Nu, 0);
    return pos;
}

void PSO::updateParticles(_real ** x, _real ** y, _real ** v) {
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

void PSO::updateParticlesWithDuConstrains(_real ** x, _real ** y, _real ** v, int test) {
    _real r1, r2;
    
    if(test == 2) {
    	printf("v = w*v + c1*r1*(y-x) + c2*r2*(global_min - x)\n");
	}
    // Update Particles
    for (int i = 0; i < S; ++i) {
        if(valid_particle[i] == 1){
    	 for (int j = 0; j < Nu; ++j) {
    	 	for (int k = 0; k < Nc; ++k){
    	 		int idx = k*Nu+j;
				if(test == 1) {
				    r1 = 0.5;
				    r2 = 0.5;
				}
				else{
				    r1 = rand_real(); //random->read();
				    r2 = rand_real(); //random->read();
				    //exporter->write_var("r1", r1);
				    //exporter->write_var("r2", r2);
                }

	            _real v_ant = v[i][idx];
	            v[i][idx] = w*v[i][idx] + c1*r1*(y[i][idx]-x[i][idx]) + c2*r2*(global_min[idx] - x[i][idx]);
	            if(test == 2) {
	            	printf("%.2f = %.2f * %.2f + %.2f*%.2f*(%.2f - %.2f) + %.2f*%.2f*(%.2f - %.2f)\n", v[i][idx], w, v_ant, c1, r1, y[i][idx], x[i][idx], c2, r2, global_min[idx], x[i][idx]);
	            }
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

int PSO::execute(_real * u_curr, int iteration, _real * last_best, _real * xref, _real * uref, _real * xss, _real * uss, _real * J, Plot * graph) {   
    number_of_active_particles = S;
    itr = iteration;
    
	int best_pos;
    w0 = 0.9; // initial weight
	wf = 0.1; // final weight
	w = w0;
	slope = (wf-w0)/maxiter;
	ini_v = max_v/10;

 	random = new RandomNumber(0);

	initializeConstrains(u_curr);

    if(controled_system->getParametrized()) {
        initializeParticles();
        for (int i = 0; i < N; ++i) {
            previous_control[i] = u_curr[i];    
        }
        //print_array("prev_u", previous_control, N, 0);

        if ((kpso == 1) && (iteration > 0)) {
            initializeLastBestKPSOParameters(last_best);
        }
        _real zero[4] = {0.0, 0.0, 0.0, 0.0};
        initializeStableZeroParameters(zero, 1);
        initializeStableZeroParameters(x_max, 2);
        initializeStableZeroParameters(x_min, 3);
    }
    else{
        initializeParticlesWithDuConstrains(u_curr, uss);

        if((iteration == 0) && (controled_system->getn_U() > 1)) {
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
    }

#ifdef GNUPLOT
    if(controled_system->getParametrized() == 0) {
        graph->plot_matrix("Initial Particles", x, S, Nu*Nc, 1);
    }
#endif

    initializeBestLocalFitness();

    // Inform System Model of PSO start to transmit reference to hardware
    controled_system->pso_init = 1;

	// ITERATIVE PROCESS
	int k = 0;  // index of iteration
	while (k < maxiter) {
        
        if(controled_system->getParametrized()) {
            evaluateFitnessAndDetectLocalBestParameters(x, xref, uref, xss, uss);
            best_pos = detectGlobalMinimum(k);
            updateParticles(x, y, v);    
        }
        else{
            evaluateFitnessAndDetectLocalBest(x, xref, uref, xss, uss);     
            best_pos = detectGlobalMinimum(k);
            updateParticlesWithDuConstrains(x, y, v, 0);    
        }

#ifdef GNUPLOT
        std::string plot_title1 = "All Particles "; 
        std::string plot_title2 = "Parameters "; 
        std::string plot_text;
        if(controled_system->getParametrized()) {
            _real ** control_signal;
            if((control_signal = alloc_matrix(S, controled_system->getN()*controled_system->getn_U())) == NULL) {throw(1);}
            for (int i = 0; i < S; ++i)
            {
                controled_system->control_from_parameters(&x[i][0], previous_control, &control_signal[i][0], controled_system->getN());    
            }

            plot_text = plot_title1 + std::to_string(k);
            graph->plot_matrix(plot_text, control_signal, S, controled_system->getN()*controled_system->getn_U(), 1);
            //graph->plot_matrix("Parameters", x, S, Nu*Nc, 3);
            plot_text = plot_title2 + std::to_string(k);
            graph->plot_matrix(plot_text, x, S, Nu*Nc, 2);
            graph->plot_vector("Cost", bestfitness, k, 3);

            free_matrix(control_signal, S);
        }
        else{
            graph->plot_vector("Best Particle", global_min, Nu*Nc, 2);
            graph->plot_matrix("All particles", x, S, Nu*Nc, 3);
        }
        usleep(100000);
#endif
        
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
//    if(controled_system->getParametrized()) {
//        controled_system->control_from_parameters(global_min, previous_control, last_best, 0);
//    }
//    else{
        for (int i = 0; i < Nu*Nc; ++i) {
            last_best[i] = global_min[i];
        }        
//    }


	*J = bestfitness[k-1];

	delete(random);
	return k;
}
