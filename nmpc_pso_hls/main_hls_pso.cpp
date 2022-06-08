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
#include "hls_system.hpp"
#include "hls_pseudorand.hpp"

#include "read_from_file.hpp"
#include "aux_functions.hpp"
#include "hls_system.hpp"

/************************** Constant Definitions *****************************/

/*
 * The following constants ...
 * 
 */

//#define SNIFFBOT_CONFIG
//#define PSO_CONFIG
const char config_file_std[] = "./config/sniffbot/project_config.txt";
const char sim_config_file_std[] = "./config/sniffbot/simulation_config_ring.txt";

#ifdef PSO_CONFIG
    #define _n_S             10
    #define _KPSO            1
    #define _stable_zero     1
    #define _maxiter         200
    #define _max_v           10
    #define _w0              0.9
    #define _wf              0.4
    #define _c1              2.0
    #define _c2              1.0
    #define _threshold       1e-2
    #define _stop_criteria   0
    const _real rand_min = -1.0;
    const _real rand_max = 1.0;
    const int rand_seed[_n_S] = {
        84680,
        577726,
        273600,
        804402,
        747952
    };
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
    const _real  _du_max[] =  {
        20, 
        20, 
        20, 
        20
    };
    const _real  _uss[] = {0.0, 0.0, 0.0, 0.0};

    #define  _Parametrization 0
    const _real  _Lambda = 5;
    const _real  _q_param = 2;
    const _real  _pmax[] =  {1, 1, 1, 1};
    const _real  _pmin[] =  {-1, -1, -1, -1};

    const unsigned short _controlled_state[] = {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
    const _real  _state_upper_limits[_Nx] =  {1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3};
    const _real  _state_lower_limits[_Nx] =  {-1e3, -1e3 -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3};
    const _real  _Q[] =  {
        50, 
        50, 
        50, 
        50, 
        50, 
        50, 
        0, 0, 0, 0, 0, 0};
    const _real  _Qf[] =  {
        200, 
        200, 
        200, 
        200, 
        200, 
        200, 
        0, 0, 0, 0, 0, 0};
    const _real  _R[] =  {
        0.0005, 
        0.0005, 
        0.0005, 
        0.0005
    };

    #define  _Rising_Time 0
    const _real _tr[] =  {10, 10, 10, 0, 0, 8, 0, 0, 0, 0, 0, 0};
    const _real  initial_state[] =  {7.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif


/***************** Macros (Inline Functions) Definitions *********************/

/*
 * Wrapper for non_linear solver
 */
int nonlinear_solver_wrapper(
    volatile float *x_curr,//[_Nx], 
    volatile float *u_curr,//[_n_U], 
    int iteration, 
    volatile float *last_best,//[_Nu*_n_U], 
    volatile float *xref,//[_Nu*_Nx], 
    volatile float *uref,//[_n_U], 
    volatile float *xss,//[_Nx],
    volatile float *uss,//[_n_U], 
    
    float *new_best,//[_Nu*_n_U],
    float *J
);

/*
 * Simulation Functions
 */

void update_xss(int iter, int xss_index);
void initialize_LastBest(float * last_best);
void initializeFilteringCoefficients(int number_of_states);
void add_disturbance(int iter);

/************************** Variable Definitions *****************************/

/*
 * The following are declared globally so they are zeroed
 */


/*
 * Simulation Variables
 */
char *config_file;
int twist_ref = 0;
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
int disturbance_size = 0;
_real ** disturbance_history;
_real ** disturbance_matrix;
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
    twist_ref = 0;
    char *config_file;
    char *sim_config_file;    
//    try {
        if (argc == 3){
            config_file = argv[1];
            sim_config_file = argv[2];
        }else{
            std::cout << "Running with standard project_config file." << std::endl;
            config_file = (char *) malloc(sizeof(config_file_std)*sizeof(config_file_std));
            sim_config_file = (char *) malloc(sizeof(sim_config_file_std)*sizeof(sim_config_file_std)); 
            std::memcpy(config_file, config_file_std, sizeof(config_file_std)*sizeof(config_file_std[0]));
            std::memcpy(sim_config_file, sim_config_file_std, sizeof(sim_config_file_std)*sizeof(sim_config_file_std[0]));
        }
        // config_file = (char *) malloc(sizeof(config_file_tmp)*sizeof(config_file_tmp));
        // std::memcpy
        std::cout << "-file : " << config_file << std::endl;
        std::string file(config_file);

        std::cout << "- Opening project config file." << std::endl;
        std::ifstream _config;
        _config.open(config_file, std::ios::in);
        if(!_config.is_open()){
            std::cout << "Project file not found.";
            return -1;
        }
        std::cout << "- Reading info from project config file." << std::endl;
        sim_type = read_string((std::ifstream *)&_config, (std::string)"sim_type");
        simulation_config_file = read_string((std::ifstream *)&_config, (std::string)"simulation_config_file");
        _config.close();
                
        std::ifstream sim_config;
//        sim_config.open(simulation_config_file, std::ios::in);
        sim_config.open(sim_config_file, std::ios::in);
        if(!sim_config.is_open()){
            std::cout << "System model config file not found.";
            return -1;
        }
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
        xref = (_real *) malloc(qtd_pontos*_Nx*sizeof(_real));
        uref = (_real *) malloc(qtd_pontos*_n_U*sizeof(_real));
        xss =  (_real *) malloc(_Nx*sizeof(_real));
        disturbance = (_real *) malloc(_n_U*sizeof(_real));
        if((state_matrix = alloc_matrix(ref_size, _Nx+1)) == NULL) {return(-1);}
        if((state_history = alloc_matrix(qtd_pontos, _Nx)) == NULL) {return -1;}
        if((control_history = alloc_matrix(qtd_pontos, _n_U)) == NULL) {return -1;}
        // if((disturbance_history = alloc_matrix(qtd_pontos, _n_U+1)) == NULL) {return -1;}
        // if((disturbance_matrix = alloc_matrix(disturbance_size, _n_U+1)) == NULL) {return -1;}

        float state_tmp;
        for (int j = 0; j < ref_size; ++j) {
            for (int i = 0; i < (_Nx+1); i++) {
                sim_config >> state_tmp;
                state_matrix[j][i] = state_tmp;
            }
        }

        // Set up initial Control and State Reference vectors
        int k = 0;
        int index;
        float sim_time = 0.0;
        for (int i = 0; i < qtd_pontos; ++i) {
            int index_xref = i*_Nx;
            int index_uref = i*_n_U;
            for (int j = 0; j < _Nx; ++j) {
                xref[index_xref+j] = state_matrix[k][j+1];
            }
            for (unsigned int j = 0; j < _n_U; j++){
                uref[index_uref+j] = u_ref_input[j];
            }
            if(sim_time >= state_matrix[k][0]) {
                if(k < (ref_size-1))
                    k++;
            }
            sim_time = sim_time + _Ts;
        }
      

//	}
//    catch(std::runtime_error &e) {
//		std::cout << "Exception: " << e.what() << std::endl;
//	}
    std::cout << "- Loaded all." << std::endl;

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
	int qtd_iter = 0;

	float * last_best;
	float * new_best;
	float * curr_state;
	float * next_state;
	float J_cost;

	last_best = (float *) malloc(_Nu*_n_U*sizeof(float));
	new_best = (float *) malloc(_Nu*_n_U*sizeof(float));
	curr_state = (float *) malloc(_Nx*sizeof(float));
	next_state = (float *) malloc(_Nx*sizeof(float));
    iterations = (float *) malloc(qtd_pontos*sizeof(float));
    cost_history = (float *) malloc(qtd_pontos*sizeof(float));

    initialize_LastBest(last_best);

    // Measure Execution  Time
    struct timespec requestStart, requestEnd;
    struct timespec requestStartCycle, requestEndCycle;
    double cycle_time, max_cycle_time = 0.0;
    clock_gettime(CLOCK_ID, &requestStart);

    float u_curr[_n_U];
    float v_curr[_n_U];
    float u_ant[_n_U];

    memset_loop<float>(u_curr, (float)0.0, _n_U);
    memset_loop<float>(v_curr, (float)0.0, _n_U);
    memset_loop<float>(u_ant, (float)0.0, _n_U);
    
    for (int i = 0; i < _Nx; i++)
        curr_state[i] = initial_state[i];

//    qtd_pontos = 1;
	while(iter < qtd_pontos - _Nh) {

        update_xss(iter, xss_index);

        clock_gettime(CLOCK_ID, &requestStartCycle);

        int index_xref = iter*_Nx;
        int index_uref = iter*_n_U;
        // qtd_iter = solver->execute(curr_state, u_curr, iter, last_best, &xref[iter*_Nx], &uref[iter], xss, uss, &J_cost);
        qtd_iter = nonlinear_solver_wrapper(curr_state, u_curr, iter, last_best, &xref[index_xref], &uref[index_uref], xss, uss, new_best, &J_cost);
        
        clock_gettime(CLOCK_ID, &requestEndCycle);
        cycle_time = ( requestEndCycle.tv_sec - requestStartCycle.tv_sec ) + ( requestEndCycle.tv_nsec - requestStartCycle.tv_nsec ) / BILLION;
        if(cycle_time > max_cycle_time) {
            max_cycle_time = cycle_time;
        }

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
        
        // ===============================================
        // Insert Disturbance
        // ===============================================
        for (int j = 0; j < _n_U; ++j) {
            disturbance[j] = u_curr[j];
        }
        /*
         * No disturbance for now
        if(disturbance_size > 0){
            //std::cout << "With disturbance" << std::endl;
            add_disturbance(iter);
        }
        */
        // output model
        current_system.one_step_prediction(next_state, curr_state, disturbance);
        
		for (int i = 0; i < _Nx; ++i) {
			state_history[iter][i] = next_state[i];
        }

        // Read from sensor
        for (int i = 0; i < _Nx; i++)
            curr_state[i] = next_state[i];


		iter++;
        
        for (int i = 0; i < _n_U; ++i) {
            u_ant[i] = u_curr[i];
        }

#ifdef PRINT_TO_TERMINAL        
        // Print Iteration Output
        std::cout << "Sim iter[" << iter << "]"<< std::endl;
        std::cout << "\t State"; print_formatted_float_array(curr_state, _Nx, 2, 6);
        std::cout << std::endl;
        std::cout << "\t Set  "; print_formatted_float_array(&xref[iter*_Nx], _Nx, 2, 6); 
        std::cout << std::endl;
        std::cout << "\t Control"; print_formatted_float_array(u_curr, _n_U, 2, 6); 
        std::cout << std::endl;
        std::cout << "\t Ji   " << J_cost; 
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
    if(state_history != NULL) free(state_history);
    if(control_history != NULL) free(control_history);
    // if(disturbance_history != NULL) free(disturbance_history);
    // if(disturbance_matrix != NULL) free(disturbance_matrix);

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

void initialize_LastBest(float * last_best){
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
            std::cout << "Perturbacao" << std::endl;
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

int nonlinear_solver_wrapper(
    volatile float *x_curr,//[_Nx], 
    volatile float *u_curr,//[_n_U], 
    int iteration, 
    volatile float *last_best,//[_Nu*_n_U], 
    volatile float *xref,//[_Nu*_Nx], 
    volatile float *uref,//[_n_U], 
    volatile float *xss,//[_Nx],
    volatile float *uss,//[_n_U], 
    
    float *new_best,//[_Nu*_n_U],
    float *J
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

    // System Core Generator
    typedef System<_real,_Nh, _Nx, _n_U, _Nu> _system_t; 
    _system_t _hw_system(
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
    );
    _system_t *_hw_system_ptr = &_hw_system;

    // Pseudo Random Core Generator
    typedef pseudoRand_gen<_real, _n_S> _randCore_t;
    _randCore_t _hw_rand_core(
        (int *)rand_seed, 
        (const float)rand_min, 
        (const float)rand_max
    );
    _randCore_t *_hw_rand_core_ptr = &_hw_rand_core;

    // Nonlinear PSO Solver
    const _real slope = (_wf-_w0)/_maxiter;
    const _real slope_init = 0;

    typedef PSO<_real, _randCore_t, _system_t,_n_S, _maxiter, _Nh, _Nx, _n_U, _Nu> T_solver;
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
        _uss
        ,
        // _controlled_state,
        // _state_upper_limits, 
        // _state_lower_limits, 
        // _Q, 
        // _Qf, 
        // _R, 
        // _Ts
        // ,
        _hw_system_ptr,
        _hw_rand_core_ptr
	);

    int iterations;
	_real my_x_curr[_Nx] ;
	_real my_u_curr[_n_U] ;
	_real my_last_best[_n_U*_Nu] ;
	_real my_xref[_Nx*_Nh];
	_real my_uref[_n_U] ;
	_real my_xss[_Nx];
	_real my_uss[_n_U] ;

	_real my_new_best[_Nu*_n_U];
	_real my_J;

//#pragma HLS bind_storage variable=local_control_guess type=FIFO impl=LUTRAM

    memcpy_loop_rolled<_real, float, _Nx>(my_x_curr, (volatile float *)x_curr );
    memcpy_loop_rolled<_real, float, _n_U>(my_u_curr, (volatile float *)u_curr );
    memcpy_loop_rolled<_real, float, _Nu*_n_U>(my_last_best, (volatile float *)last_best );
    memcpy_loop_rolled<_real, float, _Nu*_Nx>(my_xref, (volatile float *)xref );
    memcpy_loop_rolled<_real, float, _n_U>(my_uref, (volatile float *)uref );
    memcpy_loop_rolled<_real, float, _Nx>(my_xss, (volatile float *)xss );
    memcpy_loop_rolled<_real, float, _n_U>(my_uss, (volatile float *)uss );
    
    iterations = my_solver.execute(
        my_x_curr, 
        my_u_curr, 
        iteration, 
        my_last_best, 
        my_xref,
        my_uref, 
        my_xss,
        my_new_best,
        &my_J
    );

    memcpy_loop_rolled<float, _real, _Nu*_n_U>(new_best, (volatile float *)my_new_best );
    J[0] = (float) my_J;
    return iterations;
}
