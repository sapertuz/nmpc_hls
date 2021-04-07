#ifndef PROJECT_HPP
#define PROJECT_HPP

#include <string>

class Project {

private:
    std::string model_type;
    std::string solver_type;
    std::string sim_type;

    std::string model_config_file;
    std::string solver_config_file;
    std::string simulation_config_file;

    std::string output_file_name = "";

public:
	Project(std::string project_file);
    void run();
    void set_output_file_name(std::string output_file_name);
    std::string get_output_file_name();

};

#endif
