#ifndef PROJECT_HPP
#define PROJECT_HPP

#include <string>

class Project {

private:
    std::string sim_type;

    std::string simulation_config_file;

    std::string output_file_name = "";

public:
	Project(std::string project_file);
    void run();
    void set_output_file_name(std::string output_file_name);
    std::string get_output_file_name();

};

#endif
