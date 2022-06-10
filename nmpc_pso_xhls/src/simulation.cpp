#include <fstream>
#include <string.h>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "simulation.hpp"
#include "aux_functions.hpp"
#include "read_from_file.hpp"

#include "wrapper_nonlinear_solver.hpp"

Simulation::Simulation(){
    std::cout << "Wrong Simulation Call." << std::endl;
    throw("Wrong Simulation Constructor");
}

Simulation::Simulation(System * model, std::string config_file) {
	std::ifstream sim_config;
    sim_config.open(config_file, std::ios::in);

    if(!sim_config.is_open())
    	throw std::runtime_error("System model config file not found.");

    std::string temp_str;

    _real * u_ref_input;
    _real * initial_state;

	this->model = model;
	
    int Nx;
    Nx = _Nx;
    int Nc;
    Nc = _n_U;
    
    // Read Configuration File
	SimulationTime = read_real(&sim_config, "SimulationTime");	

#ifdef PRINT_TO_TERMINAL
	std::cout << "Simulation Time = " << SimulationTime << std::endl;
#endif

	qtd_pontos = (int)((_real) SimulationTime / _Ts);
	initial_state = (_real *) malloc(Nx*sizeof(_real));
	read_real_vector(&sim_config, "initial_state", initial_state, _Nx);

    u_ref_input = (_real *) malloc(Nc*sizeof(_real));
    read_real_vector(&sim_config, "u_ref", u_ref_input, _n_U);

    uss =  (_real *) malloc(Nc*sizeof(_real));
    read_real_vector(&sim_config, "uss", uss, Nc);
	
    ref_size = read_int(&sim_config, "state_ref");

	// Allocate memory for state and control references
	vector_size = qtd_pontos + _N;
    xref = (_real *) malloc(vector_size*Nx*sizeof(_real));
	uref = (_real *) malloc(vector_size*sizeof(_real));
    xss =  (_real *) malloc(Nx*sizeof(_real));
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
			this->xref[i*Nx+j] = state_matrix[k][j+1];
		}
//#ifdef DEBUG
//        print_array("xref", &xref[i*Nx], _Nx, 0);
//#endif
		this->uref[i] = u_ref_input[0];
	}

	// Initial control value
	u_curr = (_real *) malloc(Nc*sizeof(_real));
	for (int i = 0; i < Nc; ++i) {
		u_curr[i] = uss[i];
	}

    // Initial disturbance value
    disturbance = (_real *) malloc(Nc*sizeof(_real));

    // Read disturbance matrix from file
    try{
#ifdef PRINT_TO_TERMINAL  
        std::cout << "Search disturbance in simulation file." << std::endl;
#endif
        find_token(sim_config, "disturbance");
        disturbance_size = read_int(&sim_config, "disturbance");
#ifdef PRINT_TO_TERMINAL  
        std::cout << "disturbance_size: " << disturbance_size << std::endl;
#endif
        //sim_config >> disturbance_size;
        if((disturbance_matrix = alloc_matrix(disturbance_size, _n_U+1)) == NULL) {throw(1);}
        for (int i = 0; i < disturbance_size; i++) {
            for (int j = 0; j < (_n_U+1); ++j) {
                sim_config >> disturbance_matrix[i][j];
                //std::cout << "[" << i << "][" << j << "]" << disturbance_matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    catch(int e){
#ifdef PRINT_TO_TERMINAL                
        std::cout << "No disturbance set." << std::endl;
#endif
    }
    // Read Friction
    try{
        find_token(sim_config, "friction_coef");
        friction_coefficient = read_real(&sim_config, "friction_coef");
        //std::cout << "friction coefficient: " << friction_coefficient << std::endl;
    }
    catch(int e){
#ifdef PRINT_TO_TERMINAL        
        std::cout << "No friction coefficient set." << std::endl;
#endif
    }

    // Verify and initialize Filtered reference	
    if(_Rising_Time){
       alpha_ref = (_real *) malloc(Nx*sizeof(_real));
       initializeFilteringCoefficients(Nx);
    }

	// Initialize memory for output statistics data
    if((state_history = alloc_matrix(qtd_pontos, _Nx)) == NULL) {throw(1);}
    if((control_history = alloc_matrix(qtd_pontos, Nc)) == NULL) {throw(1);}
	iterations = (_real *) malloc(qtd_pontos*sizeof(_real));
	cost_history = (_real *) malloc(qtd_pontos*sizeof(_real));

    if((disturbance_history = alloc_matrix(qtd_pontos, Nc)) == NULL) {throw(1);}
    

	sim_config.close();

}

Simulation::~Simulation() {
    if(alpha_ref != NULL) free(alpha_ref);
    if(xref != NULL) free(xref);
	if(uref != NULL) free(uref);
    if(xss != NULL) free(xss);
    if(uss != NULL) free(uss);
    if(state_matrix != NULL) free_matrix(state_matrix,ref_size);
	if(state_history != NULL) free_matrix(state_history, qtd_pontos);
	if(control_history != NULL) free_matrix(control_history,qtd_pontos);
	if(iterations != NULL) free(iterations);
	if(cost_history != NULL) free(cost_history);
    if(u_curr != NULL) free(u_curr);
    if(disturbance != NULL) free(disturbance);
}

// Update Reference for Steady-State xss
void Simulation::update_xss(int iter, int xss_index) {
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

void Simulation::initializeFilteringCoefficients(int number_of_states){
    for (int i = 0; i < number_of_states; ++i) {
        alpha_ref[i] = exp(-3.0*_Ts/_tr[i]);
    }
}

void Simulation::initialize_LastBest(_real * last_best){
    for (int i = 0; i < _Nu*_n_U; ++i) {
        last_best[i] = 0.0;
    }    
}

void Simulation::add_disturbance(int iter){

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

void Simulation::execute(){

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

    // Measure Execution  Time
    struct timespec requestStart, requestEnd;
    struct timespec requestStartCycle, requestEndCycle;
    double cycle_time, max_cycle_time = 0.0;
    clock_gettime(CLOCK_REALTIME, &requestStart);

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

        clock_gettime(CLOCK_REALTIME, &requestStartCycle);

        // qtd_iter = solver->execute(curr_state, u_curr, iter, last_best, &xref[iter*_Nx], &uref[iter], xss, uss, &J_cost);
        qtd_iter = run_nonlinear_solver(curr_state, u_curr, iter, last_best, &xref[iter*_Nx], &uref[iter], xss, uss, new_best, &J_cost);


        for (int i = 0; i < _n_U; ++i) {
            u_curr[i] = new_best[i*_Nu];
        }
    
        clock_gettime(CLOCK_REALTIME, &requestEndCycle);
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

        // Update ddu_max
        model->updateAcceleration(v_curr, u_curr, u_ant);
        
        // output model
        model->one_step_prediction(next_state, curr_state, disturbance);
		model->setState(next_state);
        
		for (int i = 0; i < _Nx; ++i) {
			state_history[iter][i] = next_state[i];
        }

        // Read from sensor
        for (_uchar i = 0; i < _Nx; i++)
            curr_state[i] = next_state[i];


		iter++;

#ifdef PRINT_TO_TERMINAL        
        // Print Iteration Output
        std::cout << "Sim : " << std::right << std::setw(2) << iter << " | Solv : " << qtd_iter;
        std::cout << " | u :\t";
        for (int i = 0; i < _n_U; ++i)
        {
            // printf(" | u[%d]: %.2f\t", i, u_curr[i]);
        	std::cout << std::fixed << std::setprecision(2) << std::right << std::setw(6) << u_curr[i] << "\t";
        } 
        std::cout << "| J: " << std::fixed << std::setprecision(2) << std::right << std::setw(6) << J_cost << " | "; 
        for (int i = 0; i < _Nx; i++)
            std::cout << std::right << std::setw(8) << std::fixed << std::setprecision(2) << next_state[i] << " ";

        std::cout << "|" << std::endl;
#endif


    }

    // End Execution Time Measurement
    clock_gettime(CLOCK_REALTIME, &requestEnd);
    double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec ) / BILLION;

#ifdef PRINT_TO_TERMINAL
    std::cout << "===== Timing ====================="  << std::endl;
    std::cout << "Total Time: " << accum*1000 << " ms" << std::endl;
    std::cout << "Average Time per Cycle: " << accum/qtd_pontos*1000 << " ms" << std::endl;
    std::cout << "Maximum Time per Cycle: " << max_cycle_time*1000 << " ms" << std::endl;
    std::cout << "Sampling Time: " << _Ts*1000 << " ms" << std::endl;
    std::cout << "=================================="  << std::endl;
#endif

    std::cout << "Avg.: " << accum/qtd_pontos*1000 << " ms | " << "Max: " << max_cycle_time*1000 << " ms | Ts: " << _Ts*1000 << " ms" << std::endl;


#ifdef GNUPLOT
    delete(graph);
#endif

    if(last_best != NULL) free(last_best);
    if(new_best != NULL) free(new_best);
	if(curr_state != NULL) free(curr_state);
	if(next_state != NULL) free(next_state);
}

void export_matrix_to_file(std::ofstream& file, _real ** matrix, int rows, int columns, std::string name) {
	file << name << " = [";
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			file << std::scientific << matrix[i][j] << " ";
		}
		file << ";" << std::endl << "\t";
	}
	file << "];" << std::endl;
}

void export_vector_to_file(std::ofstream& file, const _real * vector, int rows, std::string name, int row_size) {
	file << name << " = [";
	for (int i = 0; i < rows; ++i) {
		file << std::scientific << vector[i] << " ";
        if((row_size > 0) && (((i+1) % row_size) == 0))
            file << ";" << std::endl << "\t";
	}
	file << "];" << std::endl;
}

void Simulation::save(){
	using namespace std;

	_real * upper_limits;
	_real * lower_limits;

	ofstream sim_matlab;
	sim_matlab.open("matlab/nmpc_sim.m", ios::out);

#ifdef PRINT_TO_TERMINAL
	cout << "Export Simulation data do Matlab File" << endl;
#endif

	export_matrix_to_file(sim_matlab, state_history, qtd_pontos, _Nx, "NMPC_SIM.state_history");
	export_matrix_to_file(sim_matlab, control_history, qtd_pontos, _n_U, "NMPC_SIM.control_history");

    export_matrix_to_file(sim_matlab, disturbance_history, qtd_pontos, _n_U, "NMPC_SIM.disturbance_history");


	//export_vector_to_file(sim_matlab, control_history, qtd_pontos, "NMPC_SIM.control_history", 0);
	export_vector_to_file(sim_matlab, iterations, qtd_pontos, "NMPC_SIM.iterations", 0);
	export_vector_to_file(sim_matlab, cost_history, qtd_pontos, "NMPC_SIM.cost_history", 0);
    
    export_vector_to_file(sim_matlab, xref, vector_size*_Nx, "NMPC_SIM.xref", _Nx);
	export_vector_to_file(sim_matlab, uref, vector_size, "NMPC_SIM.uref", 0);

    export_vector_to_file(sim_matlab, _state_upper_limits, _Nx, "Model.upper_limits", 0);
	export_vector_to_file(sim_matlab, _state_lower_limits, _Nx, "Model.lower_limits", 0);

	sim_matlab << "Model.Ts = " << _Ts << ";" << endl;
	sim_matlab << "Model.SimulationTime = " << SimulationTime << ";" << endl;
	sim_matlab << "Model.PredictionHorizon = " << _N << ";" << endl;
	sim_matlab << "Model.ControlHorizon = " << _Nu << ";" << endl;


    //sim_matlab << endl << "gera_grafico_final(NMPC_SIM, Model);" << endl;

	sim_matlab.close();
}

_real Simulation::compute_mse(float range_percent){
	//_real * mse = (_real *) malloc(_Nx*sizeof(_real));
    _real * mse = new _real[_Nx]; 
	_real mse_u = 0.0;
    _real mse_total = 0;

    _real accum_u = 0.0;
    _real previous_control = 0.0;

    int start_point = qtd_pontos - (range_percent*qtd_pontos/100);

	for (int i = 0; i < _Nx; ++i) {
		mse[i] = 0.0;
		for (int j = start_point; j < (qtd_pontos-1); ++j){
            if((model->is_angle_state(i)) && (state_history[j][i] > M_PI)) {
                //mse[i] += pow(((M_PI*2-state_history[j][i])-xref[j*_Nx+i]), 2);
                mse[i] += pow(normalize_angle(state_history[j][i])-xref[j*_Nx+i], 2);
            }
            else{
			    mse[i] += pow((state_history[j][i]-xref[j*_Nx+i]), 2);
            }

		}
#ifdef PRINT_TO_TERMINAL        
		std::cout << "MSE[" << i+1 << "]: " << mse[i] << std::endl;
#endif
		//mse += model->get_Q_index(i)*mse_partial/qtd_pontos;
        mse_total += mse[i];
	}
	//std::cout << "MSE_States: " << mse << std::endl;
	mse_u = 0.0;
	for (int j = start_point; j < qtd_pontos; ++j){
		for (int k = 0; k < _n_U; ++k)
		{
			mse_u += pow((control_history[j][k]), 2);

            accum_u += fabs(control_history[j][k] - previous_control);
            previous_control = control_history[j][k];
		}
		
	}
#ifdef PRINT_TO_TERMINAL    
	std::cout << "MSE_u: " << mse_u << "  ACCUM_U: " << accum_u << std::endl;
#endif
	//mse += mse_partial/qtd_pontos;
    mse_total += mse_u;

    // Export MSE file
    FILE * f_mse = fopen("mse.txt", "w");
    for (int i = 0; i < _Nx; ++i) {
        fprintf(f_mse, "%f\t", mse[i]);
    }
    fprintf(f_mse, "%f\t", mse_u);
    fprintf(f_mse, "%f\n", accum_u);
    //free(mse);
    delete[] mse;
    fclose(f_mse);

	return mse_total;
}

_real ** Simulation::getControlHistory(){
	return control_history;
}

int Simulation::getQuantidadePontos(){
	return qtd_pontos;
}
