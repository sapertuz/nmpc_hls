#include <fstream>
#include <string.h>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "aux_functions.hpp"
#include "config.hpp"
#include "wrapper_nonlinear_solver.hpp"
#include "system.hpp"
#include "current_model.hpp"

#define SimulationTime 0.3
#define u_ref 0.0
#define ref_size 2
#define qtd_pontos (int)(SimulationTime/_Ts)
#define vector_size (int)(qtd_pontos+_N)

const _real initial_state[] = {0.0, 0.0, 3.1415926536, 0.0};
const _real state_matrix[][_Nx+1] = {  
    {80.0, 0.0, 0.0, 0.0, 0.0},
    {80.0, 0.0, 0.0, 3.1415926536, 0.0}
};

//
_real * xss;
_real * iterations;
_real ** control_history;
_real * cost_history;
_real ** state_history;
_real * xref;
_real * uref;
_real * uss;

void initialize_LastBest(_real * last_best){
    for (int i = 0; i < _Nu*_n_U; ++i) {
        last_best[i] = 0.0;
    }    
}

void update_xss(int iter, int xss_index) {
    for (int k = xss_index; k < ref_size; k++) {
        if((iter*_Ts) <= state_matrix[k][0]) {
            xss_index = k;
            for (int j = 0; j < _Nx; ++j) {
                xss[j] = state_matrix[xss_index][j+1];
            }
            break;
        }
    }
}

int main(){
    // Create model
    System * model;
    model = new ModelState();
    // Variables

    xss = (_real *) malloc(_Nx*sizeof(_real));
    control_history = alloc_matrix(qtd_pontos, _n_U);
    iterations = (_real *) malloc(qtd_pontos*sizeof(_real));
	cost_history = (_real *) malloc(qtd_pontos*sizeof(_real));
    state_history = alloc_matrix(qtd_pontos, _Nx);
    uss =  (_real *) malloc(_Nc*sizeof(_real));
    xref = (_real *) malloc(vector_size*_Nx*sizeof(_real));
	uref = (_real *) malloc(vector_size*sizeof(_real));

    _real u_curr[_n_U];
    for (int i = 0; i < _n_U; ++i) {
        u_curr[i] = _uss;
    }
    // Set up initial Control and State Reference vectors
	int k = 0;
	float sim_time = 0.0;
	for (int i = 0; i < vector_size; ++i) {
		if(sim_time >= state_matrix[k][0]) {
			if(k < (ref_size-1))
				k++;
		}
		sim_time = sim_time + _Ts;
		for (int j = 0; j < _Nx; ++j) {
			xref[i*_Nx+j] = state_matrix[k][j+1];
		}
		uref[i] = 0;
	}

    int iter = 0;
    int xss_index = 0;
	int qtd_iter;

	_real * last_best;
	_real * new_best;
	_real * curr_state;
	_real * next_state;
	_real J_cost;

	last_best = (_real *) malloc(_Nu*_n_U*sizeof(_real));
	new_best = (_real *) malloc(_Nu*_n_U*sizeof(_real));
	curr_state = (_real *) malloc(_Nx*sizeof(_real));
	next_state = (_real *) malloc(_Nx*sizeof(_real));

    initialize_LastBest(last_best);

    _real v_curr[4]  = {0.0, 0.0, 0.0, 0.0};
    _real u_ant[4];

    // Initialize First read sensor
    for (_uchar i = 0; i < _Nx; i++)
    	curr_state[i] = _current_state[i];

    model->setState(curr_state);

	while(iter < qtd_pontos) {

        for (int i = 0; i < _n_U; ++i) {
            u_ant[i] = u_curr[i];
        }

        update_xss(iter, xss_index);

        qtd_iter = run_nonlinear_solver(curr_state, u_curr, iter, last_best, &xref[iter*_Nx], &uref[iter], xss, uss, new_best, &J_cost);


        for (int i = 0; i < _n_U; ++i) {
            u_curr[i] = new_best[i*_Nu];
        }
    
		// Save Data for Statistics
		for (int i = 0; i < _n_U; ++i) {
			control_history[iter][i] = u_curr[i];
			//control_history[iter*_n_U+i] = u_curr[i];
		}
		iterations[iter] = qtd_iter;
		cost_history[iter] = J_cost;
        
        // Update ddu_max
        model->updateAcceleration(v_curr, u_curr, u_ant);
        
        // output model
        model->one_step_prediction(next_state, curr_state, u_curr);
		model->setState(next_state);
        
		for (int i = 0; i < _Nx; ++i) {
			state_history[iter][i] = next_state[i];
        }

        // Read from sensor
        for (_uchar i = 0; i < _Nx; i++)
            curr_state[i] = next_state[i];


		iter++;

        // Print Iteration Output
        std::cout << "Sim iter: " << std::right << std::setw(2) << iter << " | Solv iter: " << qtd_iter;
        for (int i = 0; i < _n_U; ++i)
        {
            // printf(" | u[%d]: %.2f\t", i, u_curr[i]);
        	std::cout << " | u[" << i << "]: "<< std::fixed << std::setprecision(2) << std::right << std::setw(6) << u_curr[i] << "\t";
        } 
        std::cout << "\tJ: " << std::fixed << std::setprecision(2) << std::right << std::setw(6) << J_cost << " | "; 
        for (int i = 0; i < _Nx; i++)
            std::cout << std::right << std::setw(8) << std::fixed << std::setprecision(2) << next_state[i] << " ";

        std::cout << std::endl;
    }

}
