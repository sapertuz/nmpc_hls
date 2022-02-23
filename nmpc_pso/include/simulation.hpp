#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#define PRINT_TO_TERMINAL

#include <iostream>
#include <string>
#include <iomanip>
#include "config.h"
#include "pso.hpp"
#include "system.hpp"
#include "nonlinear_solver.hpp"

class Simulation {

protected:
    _real SimulationTime;
    int qtd_pontos;

    System * model;
    NonlinearSolver * solver;

    _real * xref; // Reference trajetory for states
    _real * uref; // Reference command
    _real * xss;  // Reference states at steady-state
    _real * uss;  // Reference command at steady-state

    int ref_size;
    _real ** state_matrix; // Reference trajectory from file

    _real * u_curr;

    _real ** control_history;
    _real ** state_history;
    _real * iterations;
    _real * cost_history;

    int vector_size;

    // Filter Parameters
    _real * alpha_ref;        // Filtered step coeficient for 95% of rising time

    _real * disturbance;
    _real ** disturbance_history;
    _real ** disturbance_matrix;
    int disturbance_size;
    _real friction_coefficient;

public:
    Simulation();
    Simulation(System * model, NonlinearSolver * solver, std::string config_file);
    virtual ~Simulation();
    virtual void execute();
    void save();
    //void plot(std::string output_file_name);
    _real compute_mse(float range_percent);

    _real ** getControlHistory();
    int getQuantidadePontos();
    void execute_verify(int qtd, _real * control_sim);

protected:
    void update_xss(int iter, int xss_index);
    void initializeFilteringCoefficients(int number_of_states);
    void reference_filter(_real * set_point, _real * xref, int horizon, int number_of_states);
    void initialize_LastBest(_real * last_best);
    void add_disturbance(int iter);
    void add_friction(int iter);
};

#endif
