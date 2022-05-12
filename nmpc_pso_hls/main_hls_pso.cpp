/***************************** Include Files *********************************/
#include <iostream>
#include <math.h>
#include <cstring>
#include <string>

#include <fstream>

#include <iomanip>

#ifdef __SYNTHESIS__
#include "hls_math.h"
typedef half _real;
#else
typedef float _real;
#endif

#include "config.hpp"
#include "hls_pso.hpp"
#include "read_from_file.hpp"
#include "aux_functions.hpp"
#include "hls_system.hpp"

/************************** Constant Definitions *****************************/

/*
 * The following constants ...
 * 
 */


#define SNIFFBOT_CONFIG
#define PSO_CONFIG

#ifdef PSO_CONFIG
    #define _n_S             10
    #define _KPSO            1
    #define _stable_zero     1
    #define _maxiter         100
    #define _max_v           10
    #define _w0              0.9
    #define _wf              0.4
    #define _c1              2.0
    #define _c2              1.0
    #define _threshold       1e-2
    #define _stop_criteria   0
#endif

#ifdef INVERTED_PENDULUM_CONFIG

    #define _Nh     3       // Prediction Horizon
    #define _Nu     3       // Control Horizon
    #define _Nx     4       // Number of States
    #define _n_U    1       // Number of Inputs

    const _real  _Ts = 0.1;
    const _real  _u_max[]  = {50.0};
    const _real  _u_min[]  = {-50.0};
    const _real  _du_max[] = {50.0};
    const _real  _uss[] = {0.0};

    #define _Parametrization 0
    const _real _Lambda = 10;
    const _real _q_param = 10;
    const _real _pmax[] = {1.0};
    const _real _pmin[] = {-1.0};

    const unsigned short _controlled_state[] = {1, 1, 1, 1};
    const _real _state_upper_limits[] = {0.5, 1e3, 1e3, 1e3} ;
    const _real _state_lower_limits[] = {-0.5, -1e3, -1e3, -1e3} ;
    const _real _Q[] = {1e3, 0.0, 1e-1, 0.0};
    const _real _Qf[] = {1e4, 0.0, 1e-0, 0.0};
    const _real _R[] = {1e-4};

    #define _Rising_Time 0
    const _real _tr[] =  {0, 0, 0, 0};
    float initial_state[] = {0.0, 0.0, 3.1415926536, 0.0};
    // float x_ss[] = {0.4, 0.3, 0.2, 0.1};
#elif defined(SNIFFBOT_CONFIG)
    #define _Nh 10
    #define _Nu 10
    #define _Nx 12
    #define _n_U 4
    const _real  _Ts = 0.05;
    const _real  _u_max[] =  {100, 100, 100, 100};
    const _real  _u_min[] =  {-100, -100, -100, -100};
    const _real  _du_max[] =  {20, 20, 20, 20};
    const _real  _uss[] = {0.0, 0.0, 0.0, 0.0};

    #define  _Parametrization 0
    const _real  _Lambda = 5;
    const _real  _q_param = 2;
    const _real  _pmax[] =  {1, 1, 1, 1};
    const _real  _pmin[] =  {-1, -1, -1, -1};

    const unsigned short _controlled_state[] = {1, 1, 1, 1};
    const _real  _state_upper_limits[_Nx] =  {1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3};
    const _real  _state_lower_limits[_Nx] =  {-1e3, -1e3 -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3};
    const _real  _Q[] =  {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
    const _real  _Qf[] =  {10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0};
    const _real  _R[] =  {0.02, 0.0, 0.0, 0.0};

    #define  _Rising_Time 0
    const _real _tr[] =  {10, 10, 10, 0, 0, 8, 0, 0, 0, 0, 0, 0};
    const _real  initial_state[] =  {7.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif

/***************** Macros (Inline Functions) Definitions *********************/

/*
 * Wrapper for non_linear solver
 */
float nonlinear_solver_wrapper(
    float x_curr[_Nx], 
    float u_curr[_n_U], 
    int iteration, 
    float last_best[_Nu*_n_U], 
    float xref[_Nu*_Nx], 
    float uref[_n_U], 
    float xss[_Nx],
    float uss[_n_U], 
    float new_best[_Nu*_n_U],
    float * J
);

/*
 * Simulation Functions
 */

void update_xss(int iter, int xss_index);
void initializeFilteringCoefficients(int number_of_states);
void initialize_LastBest(float * last_best);
void add_disturbance(int iter);

/************************** Variable Definitions *****************************/

/*
 * The following are declared globally so they are zeroed
 */


/*
 * Simulation Variables
 */

_real SimulationTime;

typedef System<float, _Nh, _Nx, _n_U, _Nu> T_system;
std::string sim_type;
std::string simulation_config_file;

int qtd_pontos;
int ref_size;
int vector_size;

float * xref; // Reference trajetory for states
float * uref; // Reference command
float * xss;  // Reference states at steady-state
float * uss;  // Reference command at steady-state

_real ** state_matrix; // Reference trajectory from file

_real * u_curr;

_real ** control_history;
_real ** state_history;
_real * iterations;
_real * cost_history;

// Filter Parameters
_real * alpha_ref;        // Filtered step coeficient for 95% of rising time

_real * disturbance;
_real ** disturbance_history;
_real ** disturbance_matrix;
int disturbance_size;
_real friction_coefficient;

/*****************************************************************************/
/**
*
* The purpose of this function is to make a simulation of the nmpc controller
* using the PSO as solver
*
* @param	argv[0] string with path to config file
* @return	0 to indicate 
*
* @note		
*
******************************************************************************/

int main(int argc, char ** argv){
    std::cout << "Start KPSO NMPC" << std::endl;
    char *config_file;
    
    try {
        std::cout << "Running with standard project_config file." << std::endl;
        char config_file[] = "./config/sniffbot/project_config.txt";

        std::string file(config_file);

        std::cout << "- Opening project config file." << std::endl;
        std::ifstream _config;
        _config.open(config_file, std::ios::in);
        if(!_config.is_open())
            throw std::runtime_error("Project file not found.");
        std::cout << "- Reading info from project config file." << std::endl;
        sim_type = read_string((std::ifstream *)&_config, (std::string)"sim_type");
        simulation_config_file = read_string((std::ifstream *)&_config, (std::string)"simulation_config_file");
        _config.close();

        std::ifstream sim_config;
        sim_config.open(simulation_config_file, std::ios::in);
        if(!sim_config.is_open())
            throw std::runtime_error("System model config file not found.");
        std::string temp_str;

        _real * u_ref_input;
        _real * initial_state;

        // Read Configuration File
        float SimulationTime = read_real(&sim_config, (std::string)"SimulationTime");
#ifdef PRINT_TO_TERMINAL
	std::cout << "Simulation Time = " << SimulationTime << std::endl;
#endif
        qtd_pontos = (int)((float) SimulationTime / _Ts);

	    initial_state = (float *) malloc(_Nx*sizeof(float));
	    read_real_vector(&sim_config, (std::string)"initial_state", initial_state, _Nx);

        u_ref_input = (_real *) malloc(_n_U*sizeof(_real));
        read_real_vector(&sim_config, (std::string)"u_ref", u_ref_input, _n_U);

        uss =  (_real *) malloc(_n_U*sizeof(_real));
        read_real_vector(&sim_config, (std::string)"uss", uss, _n_U);
        
        ref_size = read_int(&sim_config, (std::string)"state_ref");

        // Allocate memory for state and control references
        vector_size = qtd_pontos + _Nh;
        xref = (_real *) malloc(vector_size*_Nx*sizeof(_real));
        uref = (_real *) malloc(vector_size*sizeof(_real));
        xss =  (_real *) malloc(_Nx*sizeof(_real));
        if((state_matrix = alloc_matrix(ref_size, _Nx+1)) == NULL) {throw(1);}

        for (int i = 0; i < ref_size; i++) {
            for (int j = 0; j < (_Nx+1); ++j) {
    		sim_config >> state_matrix[i][j];
            }
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
    //#ifdef DEBUG
    //        print_array("xref", &xref[i*Nx], _Nx, 0);
    //#endif
            uref[i] = u_ref_input[0];
        }

	}
    catch(std::runtime_error &e) {
		std::cout << "Exception: " << e.what() << std::endl;
	}

    T_system current_system(
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
        _Ts);

    int iter = 0;
    int xss_index = 0;
	int qtd_iter;

	float * last_best;
	float * new_best;
	float * curr_state;
	float * next_state;
	float J_cost;

	last_best = (float *) malloc(_Nu*_n_U*sizeof(float));
	new_best = (float *) malloc(_Nu*_n_U*sizeof(float));
	curr_state = (float *) malloc(_Nx*sizeof(float));
	next_state = (float *) malloc(_Nx*sizeof(float));

    initialize_LastBest(last_best);

    // Measure Execution  Time
    struct timespec requestStart, requestEnd;
    struct timespec requestStartCycle, requestEndCycle;
    double cycle_time, max_cycle_time = 0.0;
    clock_gettime(CLOCK_ID, &requestStart);

    float u_curr[_n_U];
    float v_curr[_n_U];
    float u_ant[_n_U];

    memset(u_curr, 0, _n_U * sizeof(float));
    memset(v_curr, 0, _n_U * sizeof(float));
    memset(u_ant, 0, _n_U * sizeof(float));

	while(iter < qtd_pontos) {

        for (int i = 0; i < _n_U; ++i) {
            u_ant[i] = u_curr[i];
        }

        update_xss(iter, xss_index);

        clock_gettime(CLOCK_ID, &requestStartCycle);

        // qtd_iter = solver->execute(curr_state, u_curr, iter, last_best, &xref[iter*_Nx], &uref[iter], xss, uss, &J_cost);
        qtd_iter = nonlinear_solver_wrapper(curr_state, u_curr, iter, last_best, &xref[iter*_Nx], &uref[iter], xss, uss, new_best, &J_cost);


        for (int i = 0; i < _n_U; ++i) {
            u_curr[i] = new_best[i*_Nu];
        }
    
        clock_gettime(CLOCK_ID, &requestEndCycle);
        cycle_time = ( requestEndCycle.tv_sec - requestStartCycle.tv_sec ) + ( requestEndCycle.tv_nsec - requestStartCycle.tv_nsec ) / BILLION;
        if(cycle_time > max_cycle_time) {
            max_cycle_time = cycle_time;
        }

		// Save Data for Statistics
		for (int i = 0; i < _n_U; ++i) {
			control_history[iter][i] = u_curr[i];
			//control_history[iter*_n_U+i] = u_curr[i];
		}
		iterations[iter] = qtd_iter;
		cost_history[iter] = J_cost;
        
        // ===============================================
        // Insert Disturbance
        // ===============================================
        for (int j = 0; j < _n_U; ++j) {
                disturbance[j] = u_curr[j];
        }
        if(disturbance_size > 0){
            //std::cout << "With disturbance" << std::endl;
            add_disturbance(iter);
        }

        // output model
        current_system.one_step_prediction(next_state, curr_state, disturbance);
        
		for (int i = 0; i < _Nx; ++i) {
			state_history[iter][i] = next_state[i];
        }

        // Read from sensor
        for (int i = 0; i < _Nx; i++)
            curr_state[i] = next_state[i];


		iter++;

#ifdef PRINT_TO_TERMINAL        
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
#endif


    }

    // End Execution Time Measurement
    clock_gettime(CLOCK_ID, &requestEnd);
    double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec ) / BILLION;

    std::cout << "Avg.: " << accum/qtd_pontos*1000 << " ms | " << "Max: " << max_cycle_time*1000 << " ms | Ts: " << _Ts*1000 << " ms" << std::endl;

    if(last_best != NULL) free(last_best);
    if(new_best != NULL) free(new_best);
	if(curr_state != NULL) free(curr_state);
	if(next_state != NULL) free(next_state);

    return 0;
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

void initializeFilteringCoefficients(int number_of_states){
    for (int i = 0; i < number_of_states; ++i) {
        alpha_ref[i] = exp(-3.0*_Ts/_tr[i]);
    }
}

void initialize_LastBest(_real * last_best){
    for (int i = 0; i < _Nu*_n_U; ++i) {
        last_best[i] = 0.0;
    }    
}

void add_disturbance(int iter){

    float duration = 0.1; // 100 ms

    int disturbance_end;

    for (int i = 0; i < disturbance_size; ++i)
    {
        int disturbance_start = (int) disturbance_matrix[i][0] / _Ts;
        int end = (duration/_Ts) + 1; 
        disturbance_end = disturbance_start + end;
        //if((iter  > disturbance_start[i]) && (iter < disturbance_end)){
        
        if((iter  > disturbance_start) && (iter < disturbance_end)){
#ifdef PRINT_TO_TERMINAL                    
            std::cout << "Perturbação" << std::endl;
#endif            
            for (int j = 0; j < _n_U; ++j) {
                disturbance[j] = disturbance_matrix[i][j+1] + u_curr[j];
                disturbance_history[iter][j] = disturbance_matrix[i][j+1];
            }
        }
        else{
            for (int j = 0; j < _n_U; ++j) {
                disturbance[j] = u_curr[j];
            }
        }        
    }
}

float nonlinear_solver_wrapper(
    float x_curr[_Nx], 
    float u_curr[_n_U], 
    int iteration, 
    float last_best[_Nu*_n_U], 
    float xref[_Nu*_Nx], 
    float uref[_n_U], 
    float xss[_Nx],
    float uss[_n_U], 
    float new_best[_Nu*_n_U],
    float * J
){
#pragma HLS INTERFACE s_axilite port=return         bundle=control

#pragma HLS INTERFACE s_axilite port=x_curr         bundle=control
#pragma HLS INTERFACE s_axilite port=u_curr         bundle=control
#pragma HLS INTERFACE s_axilite port=iteration      bundle=control
#pragma HLS INTERFACE s_axilite port=last_best      bundle=control
#pragma HLS INTERFACE s_axilite port=xref           bundle=control
#pragma HLS INTERFACE s_axilite port=uref           bundle=control
#pragma HLS INTERFACE s_axilite port=xss            bundle=control
#pragma HLS INTERFACE s_axilite port=uss            Bundle=control
#pragma HLS INTERFACE s_axilite port=new_best       bundle=control

#pragma HLS INTERFACE m_axi depth=12    port=x_curr     offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=4     port=u_curr     offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=1     port=iteration  offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=40    port=last_best  offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=120   port=xref       offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=4     port=uref       offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=12    port=xss        offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=4     port=uss        offset=slave bundle=input
#pragma HLS INTERFACE m_axi depth=120   port=new_best   offset=slave bundle=input

    const _real slope = (_wf-_w0)/_maxiter;
    const _real slope_init = 0;

    typedef PSO<_real, _n_S, _maxiter, _Nh, _Nx, _n_U, _Nu> T_solver;
    T_solver my_solver(
        _stable_zero,
        _max_v,
        _w0,
        _wf,
        slope,
        _c1,
        _c2,
        slope_init,
        _u_min,
        _u_max,
        _du_max,
        _controlled_state,
        _state_upper_limits, 
        _state_lower_limits, 
        _Q, 
        _Qf, 
        _R, 
        _uss,
        _Ts
	);

	volatile _real my_x_curr[_Nx] ;
	volatile _real my_u_curr[_n_U] ;
	volatile _real my_last_best[_Nu*_n_U] ;
	volatile _real my_xref[_Nu*_N];
	volatile _real my_uref[_n_U] ;
	volatile _real my_xss[_Nx];
	volatile _real my_uss[_n_U] ;
	
	_real my_new_best[_Nu*_n_U];
	_real * my_J;

//#pragma HLS bind_storage variable=local_control_guess type=FIFO impl=LUTRAM

    memcpy_loop_rolled<_real, float, _Nx>(my_x_curr, x_curr );
    memcpy_loop_rolled<_real, float, _n_U>(my_u_curr, u_curr );
    memcpy_loop_rolled<_real, float, _Nu*_n_U>(my_last_best, last_best );
    memcpy_loop_rolled<_real, float, _Nu*_N>(my_xref, xref );
    memcpy_loop_rolled<_real, float, _n_U>(my_uref, uref );
    memcpy_loop_rolled<_real, float, _Nx>(my_xss, xss );
    memcpy_loop_rolled<_real, float, _n_U>(my_uss, uss );
    
    _real cf;
    my_solver.execute(
        my_x_curr, 
        my_u_curr, 
        iteration, 
        my_last_best, 
        my_xref,
        my_uref, 
        my_xss,
        my_new_best,
        my_J
    );
    return (float) cf;
}