#ifndef PLOT_HPP
#define PLOT_HPP

#include <string>
#include "config.h"
#include <float.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>

class Plot {

private:

    //FILE * file_initial;// = fopen("plot_init.dat", "w");
    //FILE * file_best; //= fopen("plot_best.dat", "w");
    //FILE * file_particles;
    FILE * gnuplotPipe;// = popen ("gnuplot -persistent", "w");
    int position;
    int max_position;

public:
	Plot(int number_of_graphs);
    ~Plot();
    void plot_matrix(std::string title, _real ** matrix, int rows, int cols, int position);
    void plot_vector(std::string title, _real * vector, int rows, int position);

private:
	void save_matrix(std::string file_name, _real ** matrix, int rows, int cols);
	void save_vector(std::string file_name, _real * vector, int rows);
	void set_position(int position);
	void update_position();
	void clear_plot(int position);
};

#endif
