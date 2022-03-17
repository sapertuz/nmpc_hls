#ifndef READ_FROM_FILE_HPP
#define READ_FROM_FILE_HPP

#include <string>
#include <fstream>
#include <iostream>

/*
int read_int(std::ifstream * file, std::string str);
int read_int_pos(std::ifstream * file, int position);
float read_real(std::ifstream * file, std::string str);
std::string read_string(std::ifstream * file, std::string str);
void read_real_vector(std::ifstream * file, std::string str, float * value, int size);
void read_int_vector(std::ifstream * file, int * value, int size, int with_comments);
int validate_string(std::string str1, std::string str2);
void read_line(std::ifstream * file, int number_of_lines);
void find_token(std::ifstream& file, std::string str);
*/

/*
#include "read_from_file.hpp"

#include <string>
#include <fstream>
*/

int validate_string(std::string str1, std::string str2) {
    if(str1.compare(str2)) {
#ifdef DEBUG
    std::cout << "---- " << str2 << " term not found (input was " << str1 << ")" << std::endl;
#endif        
        return -1;
    }
    return 0;
}

void read_line(std::ifstream * file, int number_of_lines) {
    std::string line;
    for (int i = 0; i < number_of_lines; i++) {
        std::getline(*file, line);
    }
    
}

int read_int(std::ifstream * file, std::string str) {
	std::string temp_str;
	int tmp_value;
#ifdef DEBUG
    std::cout << "-- Reading " << str << std::endl;
#endif
	*file >> temp_str >> tmp_value;
    if(validate_string(temp_str, str)) throw(1);

    return tmp_value;
}

int read_int_pos(std::ifstream * file, int position) {
    std::string temp_str;
    int tmp_value;
    
    for (int i = 0; i < position-1; i++) {
        *file >> temp_str;
    }
    *file >> tmp_value;
    return tmp_value;
}

float read_real(std::ifstream * file, std::string str) {
	std::string temp_str;
	float tmp_value;
#ifdef DEBUG
    std::cout << "-- Reading " << str << std::endl;
#endif
	*file >> temp_str >> tmp_value;
    if(validate_string(temp_str, str)){
#ifdef DEBUG
        std::cout << "---- Enter if read " << std::endl;
#endif  
        throw(1);
    } 

    return tmp_value;
}

std::string read_string(std::ifstream * file, std::string str){
    std::string temp_str;
    std::string return_str;
#ifdef DEBUG
    std::cout << "-- Reading string: " << str << std::endl;
#endif
    *file >> temp_str >> return_str;
    if(validate_string(temp_str, str)) throw(1);

    return return_str;
}

void read_real_vector(std::ifstream * file, std::string str, float * value, int size) {
	std::string temp_str;
#ifdef DEBUG
    std::cout << "-- Reading " << str << std::endl;
#endif
	*file >> temp_str;
    if(validate_string(temp_str, str)) throw(1);

    for (int i = 0; i < size; ++i) {
    	*file >> value[i];
    }
}

void read_int_vector(std::ifstream * file, int * value, int size, int with_comments){
    for (int i = 0; i < size; ++i) {
        *file >> value[i];
    }
    if(with_comments)
        read_line(file, 1);
}



void find_token(std::ifstream& file, std::string str) {
    std::string temp_str;
    int file_pos;
    int searching = 1;
    file.seekg (0, file.beg);
    while(searching) {
        file_pos = file.tellg();
        file >> temp_str;
#ifdef DEBUG        
        std::cout << "searching: " << temp_str << std::endl;
#endif
        searching = validate_string(temp_str, str);
        if(file.peek()==EOF) {
            throw 1;
        }
    }
    file.seekg (file_pos);
    
}

#endif
