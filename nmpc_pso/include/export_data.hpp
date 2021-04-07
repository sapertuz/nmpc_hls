#ifndef EXPORT_DATA_HPP
#define EXPORT_DATA_HPP

#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>

class ExportData {

private:
    FILE * data_out;

public:
	ExportData(std::string name);
    ~ExportData();
    void write_matrix(std::string name, int rows, int cols, _real ** x);
    void write_vector(std::string name, int rows, _real * x);
    void write_var(std::string name, _real x);
    void write_tag(std::string tag);
    void line_break();
};

#endif