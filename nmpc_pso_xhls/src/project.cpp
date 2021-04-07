#include <iostream>
#include "config.hpp"
#include "read_from_file.hpp"
#include "aux_functions.hpp"

#include "project.hpp"
#include "simulation.hpp"
#include "system.hpp"

Project::Project(std::string project_file){
	// Accessing config file
	std::ifstream sim_config;
    sim_config.open(project_file, std::ios::in);

    if(!sim_config.is_open())
    	throw std::runtime_error("Project file not found.");

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

    // Create model
    System * model;
    model = new ModelState();

    Simulation * sim = new Simulation(model, simulation_config_file);  

    sim->execute();
#ifndef DEBUG
        sim->save();
#endif

    _real mse = sim->compute_mse(50.0);

    delete(sim);
    delete(model);

}
