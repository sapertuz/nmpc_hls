#include <iostream>
#include <iomanip>
#include "config.h"
//#include "inverted_pendulum.hpp"
#include "pso.hpp"
#include "simulation.hpp"
#include <math.h>
#include "aux_functions.hpp"
#include "nonlinear_solver.hpp"
#include "project.hpp"
#include <string>

#include "fast_sin_cos.h"

using namespace std;

int main(int argc, char ** argv) {
	
	char * config_file;
	Project * projeto;

	initSinCosTable();

	try {
		if(argc > 2) {
			config_file = argv[1];
			string file(config_file);
			projeto = new Project(file);
			// Get output file names		
			char * output_file_name;
			output_file_name = argv[2];
			string file_out(output_file_name);						
			projeto->set_output_file_name(file_out);
		}
		else if(argc > 1) {
			config_file = argv[1];
			string file(config_file);
			projeto = new Project(file);
		}
		else {
			cout << "Running with standard project_config file." << endl;
			projeto = new Project("../config/sniffbot/project_config.txt");
		}

		projeto->run();

		//delete(projeto);

	}
	catch(runtime_error &e) {
		cout << "Exception: " << e.what() << endl;
	}

	return 0;
}
