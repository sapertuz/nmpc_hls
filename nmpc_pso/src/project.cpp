#include <iostream>
#include "config.h"
#include "read_from_file.hpp"
#include "aux_functions.hpp"

#include "project.hpp"
#include "simulation.hpp"
#include "nonlinear_solver.hpp"

#include "pso.hpp"
#include "svm.hpp"

#include "models/inverted_pendulum.hpp"
#include "models/vant.hpp"
#include "models/sniffbot.hpp"

Project::Project(std::string project_file){
	// Accessing config file
	std::ifstream sim_config;
    sim_config.open(project_file, std::ios::in);

    if(!sim_config.is_open())
    	throw std::runtime_error("Project file not found.");

    model_type = read_string(&sim_config, "model");
    model_config_file = read_string(&sim_config, "model_config_file");

    solver_type = read_string(&sim_config, "solver");
    solver_config_file = read_string(&sim_config, "solver_config_file");

    sim_type = read_string(&sim_config, "sim_type");
    simulation_config_file = read_string(&sim_config, "simulation_config_file");

    sim_config.close();
}

void Project::set_output_file_name(std::string output_file_name){
    this->output_file_name = output_file_name;
}
std::string Project::get_output_file_name(){
    return output_file_name;
}

using namespace std;
void Project::run(){

    System * model;

	if(model_type.compare("InvertedPendulum") == 0) {
		model = new InvertedPendulum(model_config_file);
        // std::cout << "Model: Inverted Pendulum" << std::endl;
    }
    else if (model_type.compare("Vant") == 0) {
        model = new Vant(model_config_file);
    }
    else if (model_type.compare("Sniffbot") == 0) {
        model = new Sniffbot(model_config_file);
    }
    else {
        throw std::runtime_error("Model type not found in Project Config File.");
    }

    NonlinearSolver * solver;
    
	if(solver_type.compare("pso") == 0) {
		solver = new PSO(model, solver_config_file);
        // std::cout << "PSO solver selected" << std::endl;
    }
    else if(solver_type.compare("svm") == 0) {
        solver = new SVM(solver_config_file, model);
        // std::cout << "SVM solver selected" << std::endl;
    }
	else {
        throw std::runtime_error("Solver type not found in Project Config File.");
	}

    if(sim_type.compare("simulation") == 0) {
        Simulation * sim = new Simulation(model, solver, simulation_config_file);  

        sim->execute();
#ifndef DEBUG
        sim->save();
#endif
#ifdef PLOT
        sim->plot("");
#endif
#ifdef SAVE_PLOT
        sim->plot(this->output_file_name);
#endif

        _real mse = sim->compute_mse(50.0);
        //std::cout << "MSE: " << mse << std::endl;

        delete(sim);
    }    
    else {
        throw std::runtime_error("Simulation type not found in Project Config File.");
    }

    // _real ** control_used;
    // _real ** control_sim;
    // int qtd;
    // qtd = sim->getQuantidadePontos();

    // if((control_sim = alloc_matrix(qtd, model->getn_U())) == NULL) {throw(1);}

    // control_used = sim->getControlHistory();

    // for (int i = 0; i < qtd; ++i) {
    //     for (int j = 0; j < model->getn_U(); ++j)
    //     {
    //         control_sim[i][j] = control_used[i][j];
    //     }
        
    // }

    delete(solver);
    delete(model);

}
