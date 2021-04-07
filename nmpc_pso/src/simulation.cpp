#include <fstream>
#include <string.h>
#include <string>
#include "simulation.hpp"
#include "aux_functions.hpp"
#include "plot.hpp"
#include "read_from_file.hpp"
#include <math.h>
#include "simulation.hpp"
#include <iostream>

Simulation::Simulation(){
    std::cout << "Wrong Simulation Call." << std::endl;
    throw("Wrong Simulation Constructor");
}

Simulation::Simulation(System * model, NonlinearSolver * solver, std::string config_file) {
	std::ifstream sim_config;
    sim_config.open(config_file, std::ios::in);

    if(!sim_config.is_open())
    	throw std::runtime_error("System model config file not found.");

    std::string temp_str;

    _real * u_ref_input;
    _real * initial_state;

	this->model = model;
	this->solver = solver;

    int Nx;
    Nx = model->getNx();
    int Nc;
    Nc = model->getn_U();
    
    // Read Configuration File
	SimulationTime = read_real(&sim_config, "SimulationTime");	

#ifdef PRINT_TO_TERMINAL
	std::cout << "Simulation Time = " << SimulationTime << std::endl;
#endif

	qtd_pontos = (int) (SimulationTime / model->getTs());
	initial_state = (_real *) malloc(Nx*sizeof(_real));
	read_real_vector(&sim_config, "initial_state", initial_state, model->getNx());
    
    u_ref_input = (_real *) malloc(Nc*sizeof(_real));
    read_real_vector(&sim_config, "u_ref", u_ref_input, model->getn_U());

    uss =  (_real *) malloc(Nc*sizeof(_real));
    read_real_vector(&sim_config, "uss", uss, Nc);
	
    ref_size = read_int(&sim_config, "state_ref");

	// Allocate memory for state and control references
	vector_size = qtd_pontos + model->getN();
    xref = (_real *) malloc(vector_size*Nx*sizeof(_real));
	uref = (_real *) malloc(vector_size*sizeof(_real));
    xss =  (_real *) malloc(Nx*sizeof(_real));
	if((state_matrix = alloc_matrix(ref_size, model->getNx()+1)) == NULL) {throw(1);}

    for (int i = 0; i < ref_size; i++) {
    	for (int j = 0; j < (model->getNx()+1); ++j) {
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
		sim_time = sim_time + model->getTs();
		for (int j = 0; j < model->getNx(); ++j) {
			this->xref[i*Nx+j] = state_matrix[k][j+1];
		}
//#ifdef DEBUG
//        print_array("xref", &xref[i*Nx], model->getNx(), 0);
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
        if((disturbance_matrix = alloc_matrix(disturbance_size, model->getn_U()+1)) == NULL) {throw(1);}
        for (int i = 0; i < disturbance_size; i++) {
            for (int j = 0; j < (model->getn_U()+1); ++j) {
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
    if(model->getFilterReference()){
       alpha_ref = (_real *) malloc(Nx*sizeof(_real));
       initializeFilteringCoefficients(Nx);
    }

	// Initialize memory for output statistics data
    if((state_history = alloc_matrix(qtd_pontos, model->getNx())) == NULL) {throw(1);}
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
        if((iter*model->getTs()) <= state_matrix[k][0]) {
            xss_index = k;
            for (int j = 0; j < model->getNx(); ++j) {
                xss[j] = state_matrix[xss_index][j+1];
            }
            break;
        }
    }
}

void Simulation::initializeFilteringCoefficients(int number_of_states){
    for (int i = 0; i < number_of_states; ++i) {
        alpha_ref[i] = exp(-3.0*model->getTs()/model->getRisingTime(i));
    }
}

// Filter for State Reference Trajectories
void Simulation::reference_filter(_real * set_point, _real * xref, int horizon, int number_of_states){
    for (int i = 0; i < horizon; ++i) {
        for (int j = 0; j < number_of_states; ++j) {
            xref[i*number_of_states+j] = set_point[j] + pow(alpha_ref[j],i)*(model->getState(j)-set_point[j]);
        }
    }
}

void Simulation::initialize_LastBest(_real * last_best){
    for (int i = 0; i < model->getNu()*model->getn_U(); ++i) {
        last_best[i] = 0.0;
    }    
}

void Simulation::add_disturbance(int iter){

    float duration = 0.1; // 100 ms

    int disturbance_end;

    for (int i = 0; i < disturbance_size; ++i)
    {
        int disturbance_start = (int) disturbance_matrix[i][0] / model->getTs();
        int end = (duration/model->getTs()) + 1; 
        disturbance_end = disturbance_start + end;
        //if((iter  > disturbance_start[i]) && (iter < disturbance_end)){
        
        if((iter  > disturbance_start) && (iter < disturbance_end)){
#ifdef PRINT_TO_TERMINAL                    
            std::cout << "Perturbação" << std::endl;
#endif            
            for (int j = 0; j < model->getn_U(); ++j) {
                disturbance[j] = disturbance_matrix[i][j+1] + u_curr[j];
                disturbance_history[iter][j] = disturbance_matrix[i][j+1];
            }
        }
        else{
            for (int j = 0; j < model->getn_U(); ++j) {
                disturbance[j] = u_curr[j];
            }
        }        
    }
}

void Simulation::add_friction(int iter){
    for (int j = 0; j < model->getn_U(); ++j) {
        disturbance[j] -= friction_coefficient * model->getState(1);
        disturbance_history[iter][j] -= friction_coefficient * model->getState(1);
    }
}

void Simulation::execute(){

	int iter = 0;
    int xss_index = 0;
	int qtd_iter;

#ifdef GNUPLOT	
	Plot * graph = new Plot(3);
#else
	Plot * graph = NULL;
#endif

	_real * last_best;
	_real * curr_state;
	_real * next_state;
	_real J_cost;

	last_best = (_real *) malloc(model->getNu()*model->getn_U()*sizeof(_real));
	curr_state = (_real *) malloc(model->getNx()*sizeof(_real));
	next_state = (_real *) malloc(model->getNx()*sizeof(_real));

    initialize_LastBest(last_best);

    // Measure Execution  Time
    struct timespec requestStart, requestEnd;
    struct timespec requestStartCycle, requestEndCycle;
    double cycle_time, max_cycle_time = 0.0;
    clock_gettime(CLOCK_REALTIME, &requestStart);

    _real v_curr[4]  = {0.0, 0.0, 0.0, 0.0};
    _real u_ant[4];

	while(iter < qtd_pontos) {

        for (int i = 0; i < model->getn_U(); ++i) {
            u_ant[i] = u_curr[i];
        }

        update_xss(iter, xss_index);

        clock_gettime(CLOCK_REALTIME, &requestStartCycle);

        if(model->getFilterReference()){
            reference_filter(xss, xref, model->getN(), model->getNx());
            qtd_iter = solver->execute(u_curr, iter, last_best, xref, &uref[iter], xss, uss, &J_cost, graph);
        }
        else{
            qtd_iter = solver->execute(u_curr, iter, last_best, &xref[iter*model->getNx()], &uref[iter], xss, uss, &J_cost, graph);   
        }

        if(model->getParametrized()){
            _real control_signal[model->getN()*model->getn_U()];
            model->control_from_parameters(last_best, u_curr, control_signal, model->getN());
            for (int i = 0; i < model->getn_U(); ++i) {
                u_curr[i] = control_signal[i*model->getN()+1];
            }
        }
        else{
            for (int i = 0; i < model->getn_U(); ++i) {
                u_curr[i] = last_best[i*model->getNu()];
            }
        }
        
        clock_gettime(CLOCK_REALTIME, &requestEndCycle);
        cycle_time = ( requestEndCycle.tv_sec - requestStartCycle.tv_sec ) + ( requestEndCycle.tv_nsec - requestStartCycle.tv_nsec ) / BILLION;
        if(cycle_time > max_cycle_time) {
            max_cycle_time = cycle_time;
        }

        // Update ddu_max
        model->updateAcceleration(v_curr, u_curr, u_ant);

		// Save Data for Statistics
		for (int i = 0; i < model->getn_U(); ++i) {
			control_history[iter][i] = u_curr[i];
			//control_history[iter*model->getn_U()+i] = u_curr[i];
		}
		iterations[iter] = qtd_iter;
		cost_history[iter] = J_cost;

		for (int i = 0; i < model->getNx(); ++i) {
			state_history[iter][i] = model->getState(i);
			curr_state[i] = model->getState(i);
		}

        // ===============================================
        // Insert Disturbance
        // ===============================================
        for (int j = 0; j < model->getn_U(); ++j) {
                disturbance[j] = u_curr[j];
        }
        if(disturbance_size > 0){
            //std::cout << "With disturbance" << std::endl;
            add_disturbance(iter);
        }
        if(friction_coefficient > 0.0){
            //std::cout << "With friction" << std::endl;
            add_friction(iter);    
        }
        
		//model->one_step_prediction(next_state, curr_state, u_curr);
        model->one_step_prediction(next_state, curr_state, disturbance);
		model->setState(next_state);

		iter++;

#ifdef PRINT_TO_TERMINAL        
        // Print Iteration Output
        std::cout << "Sim iter: " << std::right << std::setw(2) << iter << " | Solv iter: " << qtd_iter;
        for (int i = 0; i < model->getn_U(); ++i)
        {
            // printf(" | u[%d]: %.2f\t", i, u_curr[i]);
        	std::cout << " | u[" << i << "]: " << std::setprecision(2) << std::right << std::setw(6) << u_curr[i] << "\t";
        } 
        std::cout << std::setprecision(4) << std::right << std::setw(8) << "\tJ: " << J_cost << " | "; 
        for (int i = 0; i < model->getNx(); i++)
            std::cout << std::right << std::setw(8) << std::setprecision(2) << next_state[i] << " ";

        std::cout << std::endl;
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
    std::cout << "Sampling Time: " << model->getTs()*1000 << " ms" << std::endl;
    std::cout << "=================================="  << std::endl;
#endif

    std::cout << "Avg.: " << accum/qtd_pontos*1000 << " ms | " << "Max: " << max_cycle_time*1000 << " ms | Ts: " << model->getTs()*1000 << " ms" << std::endl;


#ifdef GNUPLOT
    delete(graph);
#endif

    if(last_best != NULL) free(last_best);
	if(curr_state != NULL) free(curr_state);
	if(next_state != NULL) free(next_state);
}

void Simulation::execute_verify(int qtd, _real * control_sim){
	int iter = 0;
	int qtd_iter;

	_real * last_best;
	_real * curr_state;
	_real * next_state;

	last_best = (_real *) malloc(model->getNu()*model->getn_U()*sizeof(_real));
	curr_state = (_real *) malloc(model->getNx()*sizeof(_real));
	next_state = (_real *) malloc(model->getNx()*sizeof(_real));

	for (int i = 0; i < model->getNu(); ++i) {
		last_best[i] = 0.0;
	}

	while(iter < qtd) {
        
		//qtd_iter = solver->execute(u_curr, iter, last_best, &xref[iter*model->getNx()], &uref[iter], graph);
		qtd_iter = 0;


		u_curr[0] = control_sim[iter];//last_best[0];

		// Save Data for Statistics
		control_history[iter][0] = u_curr[0];
		iterations[iter] = qtd_iter;
		for (int i = 0; i < model->getNx(); ++i) {
			state_history[iter][i] = model->getState(i);
			curr_state[i] = model->getState(i);
		}

		// Simulate system with current control value (One Step)
		model->one_step_prediction(next_state, curr_state, u_curr);

		model->setState(next_state);
        
        /*
		using namespace std;
		cout << "State: ";
		for (int i = 0; i < model->getNx(); i++) {
		   cout << next_state[i] << " ";
		}
		cout << endl << endl;
		//*/

		iter++;
        
        //std::cout << "Simulation iteration - Verify: " << iter << " | Solver iterations: " << qtd_iter << " | u: " << u_curr << std::endl;
	}

    if(last_best != NULL) free(last_best);
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

void export_vector_to_file(std::ofstream& file, _real * vector, int rows, std::string name, int row_size) {
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
	sim_matlab.open("../matlab/nmpc_sim.m", ios::out);

#ifdef PRINT_TO_TERMINAL
	cout << "Export Simulation data do Matlab File" << endl;
#endif

	export_matrix_to_file(sim_matlab, state_history, qtd_pontos, model->getNx(), "NMPC_SIM.state_history");
	export_matrix_to_file(sim_matlab, control_history, qtd_pontos, model->getn_U(), "NMPC_SIM.control_history");

    export_matrix_to_file(sim_matlab, disturbance_history, qtd_pontos, model->getn_U(), "NMPC_SIM.disturbance_history");


	//export_vector_to_file(sim_matlab, control_history, qtd_pontos, "NMPC_SIM.control_history", 0);
	export_vector_to_file(sim_matlab, iterations, qtd_pontos, "NMPC_SIM.iterations", 0);
	export_vector_to_file(sim_matlab, cost_history, qtd_pontos, "NMPC_SIM.cost_history", 0);
    
    export_vector_to_file(sim_matlab, xref, vector_size*model->getNx(), "NMPC_SIM.xref", model->getNx());
	export_vector_to_file(sim_matlab, uref, vector_size, "NMPC_SIM.uref", 0);

	upper_limits = model->getStateUpperLimits();
	lower_limits = model->getStateLowerLimits();
    export_vector_to_file(sim_matlab, upper_limits, model->getNx(), "Model.upper_limits", 0);
	export_vector_to_file(sim_matlab, lower_limits, model->getNx(), "Model.lower_limits", 0);

	sim_matlab << "Model.control_limits = [" << model->getU_min()[0] << " " << model->getU_max()[0] << " " << model->getDU_max()[0] << "];" << endl;
	sim_matlab << "Model.Ts = " << model->getTs() << ";" << endl;
	sim_matlab << "Model.SimulationTime = " << SimulationTime << ";" << endl;
	sim_matlab << "Model.PredictionHorizon = " << model->getN() << ";" << endl;
	sim_matlab << "Model.ControlHorizon = " << model->getNu() << ";" << endl;

    //sim_matlab << endl << "gera_grafico_final(NMPC_SIM, Model);" << endl;

	sim_matlab.close();

#ifdef PRINT_TO_TERMINAL
	cout << "Finished exporting" << endl;
#endif
}

/*
void Simulation::plot(std::string output_file_name){

    FILE * f_plot_r = fopen("plot_r.dat", "w");
    //FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    //FILE * gnuplotPipe = popen ("/opt/local/bin/gnuplot -persistent", "w");
    FILE * gnuplotPipe = popen ("/usr/local/opt/gnuplot/ -persistent", "w");

    std::string * plot_color;
    plot_color = new std::string[6];
    plot_color[0] = "blue";
    plot_color[1] = "dark-green";
    plot_color[2] = "black";
    plot_color[3] = "light-red";
    plot_color[4] = "orange";
    plot_color[5] = "dark-magenta";

    _real * state_upper_limits;
    _real * state_lower_limits;

    state_upper_limits = model->getStateUpperLimits();
    state_lower_limits = model->getStateLowerLimits();

    for (int i = 0; i < qtd_pontos; ++i) {
    	// Export simluation time based on sampling time
    	fprintf(f_plot_r, "%lf\t", (float) i*model->getTs());
    	// Export State History with respective reference trajectories according to plot configuration variable
    	for (int j = 0; j < model->plot_number_of_states; ++j) {
    		for (int k = 0; k < model->getNx(); ++k) {
	    		if(model->plot_states_config[j*model->getNx()+k] == 1){
	    			fprintf(f_plot_r, "%lf\t%lf\t", state_history[i][k], xref[i*model->getNx()+k]);	
	    		}
    		}
    	}
		// Export Control History with respective max and min values
    	for (int j = 0; j < model->plot_number_of_controls; ++j) {
    		for (int k = 0; k < model->getn_U(); ++k) {
	    		if(model->plot_control_config[j*model->getn_U()+k] == 1){
	    			fprintf(f_plot_r, "%lf\t%lf\t%lf\t", control_history[i][k], model->getU_max()[k], model->getU_min()[k]);
	    		}
    		}
    	}
    	fprintf(f_plot_r, "%lf", cost_history[i]);
        fprintf(f_plot_r, "\n");
    }

	fclose(f_plot_r);

    if(!output_file_name.empty()){
        fprintf(gnuplotPipe, "set term postscript enhanced color font \"Times-Roman, 10\" dashed\n set output \"%s\"\n", output_file_name.c_str());
        //gnuplot> set term postscript     (will produce postscript output)
        //gnuplot> set output "printme.ps" (output to any filename.ps you want)
        //gnuplot> replot                  (recreates plot but you don't see it, goes to file)

    }
    else{
        fprintf(gnuplotPipe, "set term %s enhanced font \"Times-Roman, 13\" dashed size 800,800\n", GNUPLOT_TERM); 
    }

    //fprintf(gnuplotPipe, "set term %s enhanced font \"Times-Roman, 13\" dashed size 800,800\n", GNUPLOT_TERM); 
    fprintf(gnuplotPipe, "set multiplot layout %d, 1 title \" %s \"\n", model->plot_number_of_states+model->plot_number_of_controls+2, model->getName().c_str()); 
    fprintf(gnuplotPipe, "set title \n set tmargin 1 \n set bmargin 0 \n set lmargin 9 \n set rmargin 2 \n unset xlabel \n unset xtics \n unset key \n");

    fprintf(gnuplotPipe, "set xrange [0:%d] \n", (int) SimulationTime);

    // Plot State History
    int state_pos = 2;
    for (int i = 0; i < model->plot_number_of_states; ++i) {
    	fprintf(gnuplotPipe, "set ylabel \"%s\" \n", model->plot_labels[i].c_str());
    	fprintf(gnuplotPipe, "plot ");
    	int color = 0;
    	for (int k = 0; k < model->getNx(); ++k) {
            if(model->plot_states_config[i*model->getNx()+k] == 1){
            	// std::cout << "Color: k: " << k << " " << plot_color[k] << std::endl;
                //fprintf(gnuplotPipe, "'plot_r.dat' using 1:%d with lines notitle linetype rgb \"blue\",", state_pos++);
                fprintf(gnuplotPipe, "'plot_r.dat' using 1:%d with lines notitle linetype rgb \"%s\",", state_pos++, plot_color[color++].c_str());
                fprintf(gnuplotPipe, "'plot_r.dat' using 1:%d with lines notitle linetype rgb \"red\" dt 2,", state_pos++);
            }
    	}
    	fprintf(gnuplotPipe, "\n");
    }
    // Plot Control History
    for (int i = 0; i < model->plot_number_of_controls; ++i) {
        fprintf(gnuplotPipe, "set title \n set xtics \n set bmargin 1 \n");
        fprintf(gnuplotPipe, "set ylabel \"%s\" \n", model->plot_labels[i+model->plot_number_of_states].c_str());
        fprintf(gnuplotPipe, "set yrange [%d:%d] \n", (int) model->getU_min()[0]-1, (int) model->getU_max()[0]+1);
        fprintf(gnuplotPipe, "set xlabel \"Time (s)\" offset 0,1 \n");
        fprintf(gnuplotPipe, "set xtics offset 0,graph 0.01 \n");
        

    	fprintf(gnuplotPipe, "plot ");
    	int color = 0;
    	for (int k = 0; k < model->getn_U(); ++k) {
            if(model->plot_control_config[i*model->getn_U()+k] == 1){
            	// std::cout << "Color: k: " << k << " " << plot_color[k] << std::endl;
                //fprintf(gnuplotPipe, "'plot_r.dat' using 1:%d with lines notitle linetype rgb \"blue\",", state_pos++);
                fprintf(gnuplotPipe, "'plot_r.dat' using 1:%d with lines notitle linetype rgb \"%s\",", state_pos++, plot_color[color++].c_str());
                fprintf(gnuplotPipe, "'plot_r.dat' using 1:%d with lines notitle linetype rgb \"red\" dt 2,", state_pos++);
                fprintf(gnuplotPipe, "'plot_r.dat' using 1:%d with lines notitle linetype rgb \"red\" dt 2,", state_pos++);
            }
    	}
    	fprintf(gnuplotPipe, "\n");
    }

    // Plot Cost Function History
    fprintf(gnuplotPipe, "set title \n set xtics \n set bmargin 0 \n");
    fprintf(gnuplotPipe, "set ylabel \"Cost\" \n");
    fprintf(gnuplotPipe, "set yrange [%d:%d] \n", 0, (int) cost_history[0]);
    fprintf(gnuplotPipe, "unset xlabel \n");
    fprintf(gnuplotPipe, "unset xtics\n");


    fprintf(gnuplotPipe, "plot ");
    fprintf(gnuplotPipe, "'plot_r.dat' using 1:%d with lines notitle linetype rgb \"%s\",", state_pos++, plot_color[0].c_str());
    fprintf(gnuplotPipe, "\n");

    // Labels
    float top_position = 0.8f;
    float left_position_0 = SimulationTime*0.01;
    float left_position_1 = SimulationTime/2 + SimulationTime*0.01;
    float left_position_2 = 3*SimulationTime/4;
    float text_offset_y = 0.12f;	

    int label_number = 1;
	fprintf(gnuplotPipe, "set label %d at %f,%f \"Model: %s\" left offset 0,.5\n", label_number++, left_position_0, top_position, model->getName().c_str());    
	fprintf(gnuplotPipe, "set label %d at %f,%f \"N: %d   Nu: %d   Ts: %.3f\" left offset 0,.5\n", label_number, left_position_0, top_position-(text_offset_y*(label_number-1)), model->getN(), model->getNu(), model->getTs());    
	for (int i = 0; i < model->getn_U(); ++i) {
		label_number++;
		fprintf(gnuplotPipe, "set label %d at %f,%f \"u-min: %.2f   u-max: %.2f   du-max: %.2f\" left offset 0,.5\n", label_number, left_position_0, top_position-(text_offset_y*(label_number-1)), model->getU_min()[i], model->getU_max()[i], model->getDU_max()[i]);    
	}

    std::string Q_str("");
    std::string Qf_str("");
    std::string R_str("");

    for (int i = 0; i < model->getNx(); ++i)
    {
        char numstr[10]; 
        sprintf(numstr, "%.0e", model->get_Q_index(i));
        Q_str = Q_str + numstr + " ";
        sprintf(numstr, "%.0e", model->get_Qf_index(i));
        Qf_str = Qf_str + numstr + " ";
    }
    for (int i = 0; i < model->getn_U(); ++i)
    {
        char numstr[10]; 
        sprintf(numstr, "%.0e", model->get_R_index(i));
        R_str = R_str + numstr + " ";
    }


    label_number++;
    fprintf(gnuplotPipe, "set label %d at %f,%f \"Q: %s \" left offset 0,.5\n", label_number, left_position_0, top_position-(text_offset_y*(label_number-1)), Q_str.c_str());    
    label_number++;
    fprintf(gnuplotPipe, "set label %d at %f,%f \"Qf: %s \" left offset 0,.5\n", label_number, left_position_0, top_position-(text_offset_y*(label_number-1)), Qf_str.c_str());    
    label_number++;
    fprintf(gnuplotPipe, "set label %d at %f,%f \"R: %s \" left offset 0,.5\n", label_number, left_position_0, top_position-(text_offset_y*(label_number-1)), R_str.c_str());    


	PSO* pso = dynamic_cast<PSO*>(solver); 
	label_number++;
	fprintf(gnuplotPipe, "set label %d at %f,%f \"PSO: \" left offset 0,.5\n", label_number++, left_position_1,top_position);    
	fprintf(gnuplotPipe, "set label %d at %f,%f \"S: %d\" left offset 0,.5\n", label_number++, left_position_1,top_position-text_offset_y, pso->getS());
	fprintf(gnuplotPipe, "set label %d at %f,%f \"KPSO: %d\" left offset 0,.5\n", label_number++, left_position_1,top_position-2*text_offset_y, pso->getKPSO());    
    fprintf(gnuplotPipe, "set label %d at %f,%f \"Stable Zero: %d\" left offset 0,.5\n", label_number++, left_position_1,top_position-3*text_offset_y, pso->getStableZero());    
	fprintf(gnuplotPipe, "set label %d at %f,%f \"Max-v: %.2f\" left offset 0,.5\n", label_number++, left_position_1,top_position-4*text_offset_y, pso->getMaxV());
	fprintf(gnuplotPipe, "set label %d at %f,%f \"Max-Iter: %d\" left offset 0,.5\n", label_number++, left_position_1,top_position-5*text_offset_y, pso->getMaxiter());
	
	fprintf(gnuplotPipe, "set label %d at %f,%f \"StopCriteria: %d\" left offset 0,.5\n", label_number++, left_position_2,top_position-text_offset_y, pso->getStopCriteria());
	fprintf(gnuplotPipe, "set label %d at %f,%f \"C1: %.1f  C2: %.1f\" left offset 0,.5\n", label_number++, left_position_2,top_position-2*text_offset_y, pso->getC1(), pso->getC2());


	//fprintf(gnuplotPipe, "unset ylabel\n unset xlabel\n set yrange [0:1] \n unset border\n unset xtics\n unset ytics\n plot 0\n");
	fprintf(gnuplotPipe, "set obj 1 rect at %f,0 size %f,2\n", SimulationTime/2, SimulationTime*0.001);
	fprintf(gnuplotPipe, "unset xlabel\n unset ylabel\n set yrange [0:1] \n unset xtics\n unset ytics\n set bmargin 1\n set tmargin 1\n plot 0\n");

    fprintf(gnuplotPipe, "unset multiplot\n");

    // Save Plot to File if output_file_name exists
    // if(!output_file_name.empty()){
    //     fprintf(gnuplotPipe, "set term postscript\n set output \"%s\"\n replot\n", output_file_name.c_str());
    //     //gnuplot> set term postscript     (will produce postscript output)
    //     //gnuplot> set output "printme.ps" (output to any filename.ps you want)
    //     //gnuplot> replot                  (recreates plot but you don't see it, goes to file)

    // }

    fprintf(gnuplotPipe, "exit \n");
    pclose(gnuplotPipe);
    delete[] plot_color;

}
*/

_real Simulation::compute_mse(float range_percent){
	//_real * mse = (_real *) malloc(model->getNx()*sizeof(_real));
    _real * mse = new _real[model->getNx()]; 
	_real mse_u = 0.0;
    _real mse_total = 0;

    _real accum_u = 0.0;
    _real previous_control = 0.0;

    int start_point = qtd_pontos - (range_percent*qtd_pontos/100);

	for (int i = 0; i < model->getNx(); ++i) {
		mse[i] = 0.0;
		for (int j = start_point; j < (qtd_pontos-1); ++j){
            if((model->is_angle_state(i)) && (state_history[j][i] > M_PI)) {
                //mse[i] += pow(((M_PI*2-state_history[j][i])-xref[j*model->getNx()+i]), 2);
                mse[i] += pow(normalize_angle(state_history[j][i])-xref[j*model->getNx()+i], 2);
            }
            else{
			    mse[i] += pow((state_history[j][i]-xref[j*model->getNx()+i]), 2);
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
		for (int k = 0; k < model->getn_U(); ++k)
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
    for (int i = 0; i < model->getNx(); ++i) {
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
