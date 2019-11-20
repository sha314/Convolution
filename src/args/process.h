//
// Created by shahnoor on 11/20/19.
//

#ifndef CONVOLUTION_PROCESS_H
#define CONVOLUTION_PROCESS_H

#include <string>
#include <vector>

bool is_number(std::string str);
//bool is_number(char* str);
std::vector<int> get_int_array(int argc, char **argv, int &index);
void test_process(int argc, char **argv);
void test_process_v1(int argc, char **argv);
void test_process_v2(int argc, char **argv);

#endif //CONVOLUTION_PROCESS_H
