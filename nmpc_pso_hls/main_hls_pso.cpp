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

#include "read_from_file.hpp"
#include "aux_functions.hpp"
#include "hls_system.hpp"

#ifdef INVERTED_PENDULUM_CONFIG
    
    #define _Nh      3      // Prediction Horizon
    #define _Nu     3      // Control Horizon
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

int qtd_pontos;
int ref_size;
int vector_size;

float * xref; // Reference trajetory for states
float * uref; // Reference command
float * xss;  // Reference states at steady-state
float * uss;  // Reference command at steady-state

float ** state_matrix; // Reference trajectory from file

typedef System<float, _Nh, _Nx, _n_U, _Nu> T_system;
std::string sim_type;
std::string simulation_config_file;

void initialize_LastBest(float * last_best){
    for (int i = 0; i < _Nu*_n_U; ++i) {
        last_best[i] = 0.0;
    }    
}

float pso_wrapper(
)
{
    return 0.0;
}

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
    clock_gettime(CLOCK_REALTIME, &requestStart);

    float u_curr[_n_U];
//    memset(u_curr, (float)0.0, _n_U*(sizeof(float)));


    // Initialize First read sensor
    // for (int i = 0; i < _Nx; i++)
    //     curr_state[i] = initial_state[i];


    //model->setState(curr_state);

    if(last_best != NULL) free(last_best);
    if(new_best != NULL) free(new_best);
	if(curr_state != NULL) free(curr_state);
	if(next_state != NULL) free(next_state);

    return 0;
}
