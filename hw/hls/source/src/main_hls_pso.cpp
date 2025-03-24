/***************************** Include Files *********************************/
#include <iostream>
#include <math.h>
#include <cstring>
#include <string>

#include <fstream>

#include <iomanip>

#include "hls_nonlinear_solver.hpp"
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
const char matlab_name_std[] = "sniffbot";

/************************** Class Definitions ****************************/


/***************** Macros (Inline Functions) Definitions *********************/

/*
 * Wrapper for non_linear solver
 */

/*
 * Simulation Functions
 */

void update_xss(int iter, int xss_index);
void initialize_LastBest(float * last_best);
void initializeFilteringCoefficients(int number_of_states);
void add_disturbance(int iter);
void update_system();

void export_matrix_to_file(
    std::ofstream& file, float ** matrix, 
    int rows, int columns, std::string name);
void export_vector_to_file(
    std::ofstream& file, const float * vector, 
    int rows, std::string name, int row_size);
void save(char * local_matlab_name);

/************************** Variable Definitions *****************************/

/*
 * The following are declared globally so they are zeroed
 */


/*
 * Simulation Variables
 */
char *config_file;
int twist_ref = 0;
float SimulationTime;

std::string sim_type;
std::string simulation_config_file;
std::string matlab_file_name;

int qtd_pontos;
int ref_size;
int vector_size;

float * xref; // Reference trajetory for states
float * uref; // Reference command
float * xss;  // Reference states at steady-state
float * uss;  // Reference command at steady-state

float ** state_matrix; // Reference trajectory from file

float * u_curr;

float ** control_history;
float ** state_history;
float ** xref_history;
float ** uref_history;
float * iterations;
float * cost_history;

// Filter Parameters
float * alpha_ref;        // Filtered step coeficient for 95% of rising time

float * disturbance;
int disturbance_size = 0;
float ** disturbance_history;
float ** disturbance_matrix;
float friction_coefficient;

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
    srand (1654785693);

    std::cout << "Start KPSO NMPC" << std::endl;
    twist_ref = 0;
    char *config_file;
    char *sim_config_file;   
    char *matlab_name; 
//    try {
        if (argc == 4){
            config_file = argv[1];
            sim_config_file = argv[2];
            matlab_name = argv[3];
        }else{
            std::cout << "Running with standard project_config file." << std::endl;
            config_file = (char *) malloc(sizeof(config_file_std)*sizeof(config_file_std));
            sim_config_file = (char *) malloc(sizeof(sim_config_file_std)*sizeof(sim_config_file_std)); 
            matlab_name = (char *) malloc(sizeof(matlab_name_std)*sizeof(matlab_name_std));
            std::memcpy(config_file, config_file_std, sizeof(config_file_std)*sizeof(config_file_std[0]));
            std::memcpy(sim_config_file, sim_config_file_std, sizeof(sim_config_file_std)*sizeof(sim_config_file_std[0]));
            std::memcpy(matlab_name, matlab_name_std, sizeof(matlab_name_std)*sizeof(matlab_name_std[0]));
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

        float * u_ref_input;
        float * initial_state;

        // Read Configuration File
        SimulationTime = read_real(&sim_config, (std::string)"SimulationTime");
#ifdef PRINT_TO_TERMINAL
	    std::cout << "Simulation Time = " << SimulationTime << std::endl;
#endif
        qtd_pontos = (int)((float) SimulationTime / _Ts);

	    initial_state = (float *) malloc(_Nx*sizeof(float));
	    read_real_vector(&sim_config, (std::string)"initial_state", initial_state, _Nx);

        u_ref_input = (float *) malloc(_n_U*sizeof(float));
        read_real_vector(&sim_config, (std::string)"u_ref", u_ref_input, _n_U);

        uss =  (float *) malloc(_n_U*sizeof(float));
        read_real_vector(&sim_config, (std::string)"uss", uss, _n_U);
        
        ref_size = read_int(&sim_config, (std::string)"state_ref");

        // Allocate memory for state and control references
        xref = (float *) malloc(qtd_pontos*_Nx*sizeof(float));
        uref = (float *) malloc(qtd_pontos*_n_U*sizeof(float));
        xss =  (float *) malloc(_Nx*sizeof(float));
        disturbance = (float *) malloc(_n_U*sizeof(float));
        if((state_matrix = alloc_matrix(ref_size, _Nx+1)) == NULL) {return(-1);}
        if((state_history = alloc_matrix(qtd_pontos, _Nx)) == NULL) {return -1;}
        if((control_history = alloc_matrix(qtd_pontos, _n_U)) == NULL) {return -1;}
        if((xref_history = alloc_matrix(qtd_pontos, _Nx)) == NULL) {return -1;}
        if((uref_history = alloc_matrix(qtd_pontos, _n_U)) == NULL) {return -1;}
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
                xref_history[i][j] = state_matrix[k][j+1];
            }
            for (unsigned int j = 0; j < _n_U; j++){
                uref[index_uref+j] = u_ref_input[j];
                uref_history[i][j] = u_ref_input[j];
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

    // top_model_t my_sim_model;
    // top_model_t *my_sim_model_ptr = &my_sim_model;

    typedef System<float, float, _Nh, _Nx, _n_U, _Nu> T_sim_system;
    T_sim_system current_system(
        (float *)_u_max,
        (float *)_u_min, 
        (float *)_du_max,
        _controlled_state,
        (float *)_state_upper_limits, 
        (float *)_state_lower_limits, 
        (float *)_Q, 
        (float *)_Qf, 
        (float *)_R, 
        (float *)_uss,
        (float)_Ts
        );

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
        // qtd_iter = nonlinear_solver_wrapper(curr_state, u_curr, iter, last_best, &xref[index_xref], &uref[index_uref], xss, uss, new_best, &J_cost);
        qtd_iter = nonlinear_solver_wrapper(
            curr_state, 
            u_curr, 
            iter, 
            last_best, 
            &xref[index_xref], 
            new_best, 
            &J_cost
        );
        
        clock_gettime(CLOCK_ID, &requestEndCycle);
        cycle_time = ( requestEndCycle.tv_sec - requestStartCycle.tv_sec ) + ( requestEndCycle.tv_nsec - requestStartCycle.tv_nsec ) / BILLION;
        if(cycle_time > max_cycle_time) {
            max_cycle_time = cycle_time;
        }

        for (int i = 0; i < _n_U; ++i) {
            // u_curr[i] = new_best[i*_Nu];
            u_curr[i] = new_best[i];
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
        
        // No disturbance for now
        // if(disturbance_size > 0){
        //     //std::cout << "With disturbance" << std::endl;
        //     add_disturbance(iter);
        // }
        
        // output model
        current_system.one_step_prediction(next_state, curr_state, disturbance);
        
		for (int i = 0; i < _Nx; ++i) {
			state_history[iter][i] = next_state[i];
        }

        // Read from sensor
        for (int i = 0; i < _Nx; i++)
            curr_state[i] = next_state[i];

        // ===============================================
        // Save New Best
        // ===============================================
        for (int j = 0; j < _Nu*_n_U; ++j) {
            last_best[j] = new_best[j];
        }
        
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
        std::cout << "\t Iter   " << qtd_iter; 
        std::cout << std::endl;
#endif

    }

    // End Execution Time Measurement
    clock_gettime(CLOCK_ID, &requestEnd);
    double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec ) / BILLION;

    std::cout << "Avg.: " << accum/qtd_pontos*1000 << " ms | " << "Max: " << max_cycle_time*1000 << " ms | Ts: " << (float)_Ts*1000 << " ms" << std::endl;

    save(matlab_name);

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

void export_matrix_to_file(
    std::ofstream& file, float ** matrix, 
    int rows, int columns, std::string name
) {
	file << name << " = [";
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			file << std::scientific << matrix[i][j] << " ";
		}
		file << ";" << std::endl << "\t";
	}
	file << "];" << std::endl;
}

void export_vector_to_file(
    std::ofstream& file, const float * vector, 
    int rows, std::string name, int row_size
) {
	file << name << " = [";
	for (int i = 0; i < rows; ++i) {
		file << std::scientific << vector[i] << " ";
        if((row_size > 0) && (((i+1) % row_size) == 0))
            file << ";" << std::endl << "\t";
	}
	file << "];" << std::endl;
}

void save(char * local_matlab_name){
	using namespace std;

	float * upper_limits;
	float * lower_limits;

	ofstream sim_matlab;
    std::string of_name = "./matlab/";
    of_name += local_matlab_name;
    of_name += ".m";
	sim_matlab.open(of_name, ios::out);

#ifdef PRINT_TO_TERMINAL
	cout << "Export Simulation data do Matlab File: " << of_name << endl;
#endif

	export_matrix_to_file(sim_matlab, state_history, qtd_pontos, _Nx, "NMPC_SIM.state_history");
	export_matrix_to_file(sim_matlab, control_history, qtd_pontos, _n_U, "NMPC_SIM.control_history");

    // export_matrix_to_file(sim_matlab, disturbance_history, qtd_pontos, _n_U, "NMPC_SIM.disturbance_history");

	// export_vector_to_file(sim_matlab, control_history, qtd_pontos, "NMPC_SIM.control_history", 0);
	export_vector_to_file(sim_matlab, iterations, qtd_pontos, "NMPC_SIM.iterations", 0);
	export_vector_to_file(sim_matlab, cost_history, qtd_pontos, "NMPC_SIM.cost_history", 0);
    
    export_matrix_to_file(sim_matlab, xref_history, qtd_pontos, _Nx, "NMPC_SIM.xref");
	export_matrix_to_file(sim_matlab, uref_history, qtd_pontos, _n_U, "NMPC_SIM.uref");

    export_vector_to_file(sim_matlab, _state_upper_limits, _Nx, "Model.x_max", 0);
	export_vector_to_file(sim_matlab, _state_lower_limits, _Nx, "Model.x_min", 0);

	export_vector_to_file(sim_matlab, _u_max, _n_U, "Model.u_max", 0);
	export_vector_to_file(sim_matlab, _u_min, _n_U, "Model.u_min", 0);

	sim_matlab << "Model.Ts = " << _Ts << ";" << endl;
	sim_matlab << "Model.SimulationTime = " << SimulationTime << ";" << endl;
	sim_matlab << "Model.PredictionHorizon = " << _Nh << ";" << endl;
	sim_matlab << "Model.ControlHorizon = " << _Nu << ";" << endl;

	sim_matlab.close();
}