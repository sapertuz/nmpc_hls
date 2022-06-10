#include <iostream>
#include <math.h>
#include <string>

#include "config.hpp"
#include "aux_functions.hpp"
#include "project.hpp"

using namespace std;

int main(int argc, char ** argv) {
	srand (1654785693);

	char * config_file;
	Project * projeto;

	cout << "Start" << endl;			
	try {
		if(argc > 1) {
			config_file = argv[1];
			string file(config_file);
			projeto = new Project(file);
		}
		else {
			cout << "Running with standard project_config file." << endl;
			projeto = new Project("../project_config.txt");
		}

		projeto->run();

		//delete(projeto);

	}
	catch(runtime_error &e) {
		cout << "Exception: " << e.what() << endl;
	}

	return 0;
}
