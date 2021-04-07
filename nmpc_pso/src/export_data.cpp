#include "export_data.hpp"

ExportData::ExportData(std::string name) {
	data_out = fopen(name.c_str(), "w");
    if (data_out == NULL) {
	    throw(1);
	}
}

ExportData::~ExportData() {
	fclose(data_out);
}

void ExportData::write_matrix(std::string name, int rows, int cols, _real ** x) {
	fprintf(data_out, "%s %d %d ", name.c_str(), rows, cols);
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			fprintf(data_out, "%f ", x[i][j]);
		}
	}
	fprintf(data_out, "\t");
}

void ExportData::write_vector(std::string name, int rows, _real * x) {
	fprintf(data_out, "%s 1 %d ", name.c_str(), rows);
	for (int i = 0; i < rows; ++i)
	{
		fprintf(data_out, "%f ", x[i]);
	}
	fprintf(data_out, "\t");
}

void ExportData::write_var(std::string name, _real x) {
	//fprintf(data_out, "%s %f\t", name.c_str(), x);
	fprintf(data_out, " %f\t", x);
}

void ExportData::write_tag(std::string tag) {
	fprintf(data_out, "%s\t", tag.c_str());	
}

void ExportData::line_break(){
	fprintf(data_out, "\n");	
}