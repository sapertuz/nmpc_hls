#include "plot.hpp"

Plot::Plot(int number_of_graphs){
    gnuplotPipe = popen ("gnuplot -persistent", "w");
    // Configure Plot
    fprintf(gnuplotPipe, "clear \n set multiplot layout %d, 1 \n", number_of_graphs);
    position = 1;
    max_position = number_of_graphs;
}

Plot::~Plot(){
	fprintf(gnuplotPipe, "unset multiplot\n exit \n");
	pclose(gnuplotPipe);
}

void Plot::save_matrix(std::string file_name, _real ** matrix, int rows, int cols) {
    FILE * file_matrix = fopen(file_name.c_str(), "w");
	for (int j = 0; j < cols; ++j) {
        fprintf(file_matrix, "%d\t", j);
        for (int i = 0; i < rows; ++i) {
            fprintf(file_matrix, "%lf\t", matrix[i][j]);
        }
        fprintf(file_matrix, "\n");
    }
    fflush(file_matrix);
    fclose(file_matrix);
}

void Plot::save_vector(std::string file_name, _real * vector, int rows) {
    FILE * file_vector = fopen(file_name.c_str(), "w");
    for (int i = 0; i < rows; ++i) {
        fprintf(file_vector, "%lf\t%lf\n", (float) i, vector[i]);
    }
    fclose(file_vector);
}

void Plot::set_position(int position) {
    if(this->position != position) {
    	while (this->position != position) {
            update_position();
            fprintf(gnuplotPipe, "set multiplot next \n");      
    	}
    }
}

void Plot::update_position() {
	if (this->position < max_position) {
		this->position += 1; 
	}
	else
		this->position = 1;	
}

void Plot::clear_plot(int position) {
    set_position(position);
    fprintf(gnuplotPipe, "clear\n");
    update_position();
    set_position(position);
}

void Plot::plot_matrix(std::string title, _real ** matrix, int rows, int cols, int position){
    std::string file_name;
    file_name = title + ".dat";

	save_matrix(file_name, matrix, rows, cols);
    clear_plot(position);

    //fprintf(gnuplotPipe, "set title \"%s\" \n unset key \n set xrange [0:%d] \n set yrange [-3:3] \n", title.c_str(), cols);
    fprintf(gnuplotPipe, "set title \"%s\" \n unset key \n set xrange [0:%d] \n", title.c_str(), cols);
    fprintf(gnuplotPipe, "plot ");
    for (int i = 0; i < rows; i++) {
        fprintf(gnuplotPipe, "'%s' using 1:%d with lines notitle, ", file_name.c_str(), i+2);
    }
    fprintf(gnuplotPipe, " \n");
    fflush(gnuplotPipe);
    update_position();
}

void Plot::plot_vector(std::string title, _real * vector, int rows, int position){
    std::string file_name;
    file_name = title + ".dat";

	save_vector(file_name, vector, rows);
    clear_plot(position);

	fprintf(gnuplotPipe, "set title \"%s\" \n unset key \nset xrange [0:%d]\n plot '%s' using 1:2 with lines notitle\n", title.c_str(), rows, file_name.c_str());
	fflush(gnuplotPipe);
	update_position();
}
