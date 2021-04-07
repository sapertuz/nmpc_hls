#ifndef SVM_HPP
#define SVM_HPP

#include <iostream>
#include <string>
#include "config.h"
#include "nonlinear_solver.hpp"
#include "export_data.hpp"

class SVM : public NonlinearSolver {

	System * controled_system;

    _real gamma;
    _real rho;
	_real * weights;
    _real * normalization;
    _real * centers;

	int number_of_inputs;
	int number_of_SVs;

	int number_of_states;
	_real * state_history;
	int * state_history_config;
	int history_size;
	int history_type;

	_real * future_references;
	int * future_references_config;
	int prediction_size;
	int prediction_type;

	_real * input_data;

	std::string training_file;

public:
	SVM(std::string config_file, System * controled_system);
	~SVM();
	int execute(_real * u_curr, int iteration, _real * last_best, _real * xref, _real * uref, _real * xss, _real * uss, _real * J, Plot * graph);
	_real calculate();

	int get_number_of_inputs();
	int get_number_of_SVs();
	int get_number_of_states();
	int get_history_size();
	int get_history_type();
	int get_predicion_size();
	int get_prediction_type();
	int * get_state_history_config();
	int * get_future_references_config();

private:
	void loadInputConfigFile(std::string config_file);
	void loadSVMTrainingData(std::string training_file);
	void allocInputVector();
	void initializeStateHistory(_real * xref);
    void updateStateHistory(_real * xref);

	void generateInputData(_real * xref);
	void normalizeInputs();
};

#endif
