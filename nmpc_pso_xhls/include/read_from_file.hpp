#ifndef READ_FROM_FILE_HPP
#define READ_FROM_FILE_HPP

#include <string>
#include <fstream>
#include "config.hpp"

int read_int(std::ifstream * file, std::string str);
int read_int_pos(std::ifstream * file, int position);
_real read_real(std::ifstream * file, std::string str);
std::string read_string(std::ifstream * file, std::string str);
void read_real_vector(std::ifstream * file, std::string str, _real * value, int size);
void read_int_vector(std::ifstream * file, int * value, int size, int with_comments);
int validate_string(std::string str1, std::string str2);
void read_line(std::ifstream * file, int number_of_lines);
void find_token(std::ifstream& file, std::string str);

#endif
