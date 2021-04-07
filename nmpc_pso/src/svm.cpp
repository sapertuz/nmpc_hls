#include <sstream>
#include <fstream>
#include <iomanip>
#include "svm.hpp"
#include "read_from_file.hpp"
#include "config.h"
#include "aux_functions.hpp"

SVM::SVM(std::string config_file, System * controled_system) {

	this->controled_system = controled_system;

    normalization = nullptr;
    
	loadInputConfigFile(config_file);
	loadSVMTrainingData(training_file);
	allocInputVector();
}

SVM::~SVM(){
	free(weights);
	free(normalization);
	free(centers);
	free(state_history);
	free(future_references);
	free(state_history_config);
	free(future_references_config);
}

int SVM::get_number_of_inputs(){return number_of_inputs;}
int SVM::get_number_of_SVs(){return number_of_SVs;}
int SVM::get_number_of_states(){return number_of_states;}
int SVM::get_history_size(){return history_size;}
int SVM::get_history_type(){return history_type;}
int SVM::get_predicion_size(){return prediction_size;}
int SVM::get_prediction_type(){return prediction_type;}
int * SVM::get_state_history_config(){return state_history_config;}
int * SVM::get_future_references_config(){return future_references_config;}

void SVM::loadInputConfigFile(std::string config_file){
	std::ifstream sim_config;
	std::string str_temp;
    sim_config.open(config_file, std::ios::in);

    if(!sim_config.is_open()) {
    	throw std::runtime_error("SVM config file not found.");
    }

	sim_config >> training_file;
    number_of_states = read_int(&sim_config, "number_of_states");
    history_size = read_int(&sim_config, "history_size");
    history_type = read_int(&sim_config, "history_type");
    read_line(&sim_config, 1);
    prediction_size = read_int(&sim_config, "prediction_size");
    prediction_type = read_int(&sim_config, "prediction_type");
    read_line(&sim_config, 3);

    state_history_config = (int *) malloc(number_of_states*history_size*sizeof(int));
    for (int i = 0; i < number_of_states; ++i) {
    	read_int_vector(&sim_config, &state_history_config[i*history_size], history_size, 1);
	}
	
	read_line(&sim_config, 2);
	future_references_config = (int *) malloc(number_of_states*prediction_size*sizeof(int));	
    for (int i = 0; i < number_of_states; ++i) {
    	read_int_vector(&sim_config, &future_references_config[i*prediction_size], prediction_size, 1);
	}

	sim_config.close();
}

void SVM::loadSVMTrainingData(std::string training_file) {
	// Accessing config file
	std::ifstream sim_config;
    sim_config.open(training_file, std::ios::in);

    if(!sim_config.is_open()) {
    	throw std::runtime_error("SVM Training Data file not found.");
    }

	std::string temp_str;
    int * indexes;
    int j;
    
    read_line(&sim_config, 6);
    gamma = read_real(&sim_config, "Gama:");
    
    read_line(&sim_config, 5);
    
    number_of_SVs = read_int_pos(&sim_config, 3);
    
    rho = read_real(&sim_config, "RHO:");

    read_line(&sim_config, 2);
    
    //std::cout << "Weights" << std::endl;
    // Read Weights and Indexes
    weights = (_real *) malloc(number_of_SVs*sizeof(_real));
    indexes = (int *) malloc(number_of_SVs*sizeof(int));
    for (int i = 0; i < number_of_SVs; i++) {
        sim_config >> indexes[i] >> weights[i];
        //std::cout << weights[i] << std::endl;
    }
    
    // Read normalization Data
    number_of_inputs = 0;
    
    find_token(sim_config, "Max");
    read_line(&sim_config, 1);
    
    long long pos = sim_config.tellg();
    
    int end_of_min = 1;
    while(end_of_min) {
        char *line = new char [200];
        sim_config.getline(line, 200);
        std::istringstream iss(line);
        if(!std::string(line).empty()) {
            number_of_inputs++;
        }
        else {
            end_of_min = 0;
        }
    }
    sim_config.seekg(pos);

    normalization = (_real *) malloc(number_of_inputs*2*sizeof(_real));
    
    for (int i = 0; i < number_of_inputs*2; i++) {
        sim_config >> normalization[i];
    }
    
    // Read Center Values
    find_token(sim_config, "treinamento");
    read_line(&sim_config, 1);
    
    centers = (_real *) malloc(number_of_SVs*number_of_inputs*sizeof(_real));
    j = 0;
    int m = 0;
    _real tmp;
    
    for(int i = 0; i < number_of_SVs; i++) {
    	//std::cout << "Index[" << i+1 << "]: " << indexes[i] << " m[" << m+1 << "]";
        if(indexes[i] == ++m) {
        	//std::cout << " <= " << std::endl;
            for (int k = 0; k < number_of_inputs; k++) {
                sim_config >> centers[j++];
                //std::cout << "[" << j-1  << "] " << centers[j-1];
            }
            sim_config >> tmp;
            //std::cout << " (" << tmp << ")" << std::endl;
        }
        else {
        	//std::cout << " not " << std::endl;
            for (int k = 0; k < (number_of_inputs + 1); k++) {
                sim_config >> tmp;
                //std::cout << " (" << tmp << ")";
            }
            //std::cout << std::endl;
            i--;
        }
    }
    /*
    std::cout << "Centers: [" << j-1 << "]" << std::endl;
    for (int i = 0; i < number_of_SVs*number_of_inputs; ++i) {
    	//for (int j = 0; j < number_of_inputs; ++j)
    	//{
    		std::cout << centers[i] << " ";
    	//}
    	if(((i+1) % number_of_inputs) == 0)
    		std::cout << std::endl;
    } //*/

    sim_config.close();
}


int SVM::execute(_real * u_curr, int iteration, _real * last_best, _real * xref, _real * uref, _real * xss, _real * uss, _real * J, Plot * graph){

	if(iteration == 0) {
		initializeStateHistory(xref);
	}

    updateStateHistory(xref);
	generateInputData(xref);

	normalizeInputs();

    last_best[0] = calculate();

    return 0;
}

void SVM::allocInputVector() {
	int size = 0;

	for (int i = 0; i < history_size*number_of_states; ++i) {
		if(state_history_config[i] == 1) {
			size++;
		}
	}

	for (int i = 0; i < prediction_size*number_of_states; ++i) {
		if(future_references_config[i] == 1) {
			size++;
		}
	}
	input_data = (_real *) malloc(size*sizeof(_real));
}

void SVM::initializeStateHistory(_real * xref){

	state_history = (_real *) malloc(history_size*number_of_states*sizeof(_real));
	future_references = (_real *) malloc(prediction_size*number_of_states*sizeof(_real));

	for (int i = 0; i < number_of_states; ++i) {
		for (int j = 0; j < history_size; ++j) {
			if(history_type) {
				state_history[i*history_size+j] = controled_system->getState(i) - xref[i];
			}
			else {
				state_history[i*history_size+j] = controled_system->getState(i);
			}
		}
	}
}

void SVM::updateStateHistory(_real * xref){
    _real * state_history_tmp;
    state_history_tmp = (_real *) malloc(history_size*number_of_states*sizeof(_real));
    
    for (int i = 0; i < number_of_states; i++) {
    	if(history_type) {
        	state_history_tmp[i*history_size] = controled_system->getState(i) - xref[i];
    	}
    	else {
    		state_history_tmp[i*history_size] = controled_system->getState(i);	
    	}
    }
    
    for (int i = 0; i < number_of_states; ++i) {
        for (int j = 1; j < history_size; ++j) {
            state_history_tmp[i*history_size+j] = state_history[i*history_size+j-1];
        }
    }
    
    delete(state_history);
    state_history = state_history_tmp;
}

void SVM::generateInputData(_real * xref){
	int k = 0;
   
    for (int i = 0; i < (number_of_states*history_size); ++i) {
        if(state_history_config[i] == 1) {
            input_data[k++] = state_history[i];
        }
    }
    
    for (int i = 0; i < number_of_states; ++i) {
        for (int j = 0; j < prediction_size; ++j) {
            if(future_references_config[i*prediction_size+j] == 1) {
            	if(prediction_type == 1) {
            		input_data[k++] = controled_system->getState(i) - xref[j*number_of_states+i];
            	}
            	else {
            		input_data[k++] = xref[j*number_of_states+i];	
            	}
                
            }
        }
    }
}

void SVM::normalizeInputs() {
    for (int i = 0; i < number_of_inputs; i++) {
        if((normalization[i*2] < 0) || (normalization[i*2+1] > 1)) {
            input_data[i] = (input_data[i] - normalization[i*2])/(normalization[i*2+1] - normalization[i*2]);
        }
    }
}

_real SVM::calculate() {
	_real phi;
	_real sv_out = 0;
    
    //std::cout << "Calcula SVM" << std::endl;
	for (int i = 0; i < number_of_SVs; ++i) {
	 	phi = 0;
	 	for (int j = 0; j < number_of_inputs; ++j) {
	 		phi = phi + pow(input_data[j]-centers[i*number_of_inputs+j],2);
            
            //std::cout << "SV[" << i << "] " << input_data[j] << " - " << centers[i*number_of_inputs+j] << " [" << i*number_of_inputs+j << "] | " << std::endl;
	 	}
	 	phi = exp(-phi*gamma);

        // std::cout << "Phi: " << phi << std::endl;

	 	sv_out += weights[i]*phi;
	}

    // std::cout << "Out: " << sv_out;

    sv_out -= rho;

    // std::cout << " SV_out: " << sv_out << std::endl;

    //std::cout << "Gamma: " << gamma << " Rho: " << rho << std::endl;
    
    return sv_out;
}
