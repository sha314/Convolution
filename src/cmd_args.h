//
// Created by shahnoor on 11/18/19.
//

#ifndef CONVOLUTION_CMD_ARGS_H
#define CONVOLUTION_CMD_ARGS_H


#include <vector>
#include <string>

namespace
{
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace


void get_option_a(int argc, char *const *argv, std::vector<int> &a_usecols, std::vector<std::string> &a_names, int i);

void get_option_b(int argc, char *const *argv, std::vector<int> &b_usecols, std::vector<std::string> &b_names, int i);

int parse_cmd_arg_boost(int argc, char *const *argv, std::string &in_filename, std::string &out_filename, std::vector<int> &a_usecols,
                         std::vector<int> &b_usecols, std::string &info, bool &write_header_and_comment, int &skiprows,
                         bool &write_input_data, int &f_precision, int &n_threads, double &threshold, int &times, char& delimiter);

void parse_cmd_arg(int argc, char *const *argv, std::string &in_filename, std::string &out_filename, std::vector<int> &a_usecols,
                   std::vector<int> &b_usecols, std::string &info, bool &write_header_and_comment, int &skiprows,
                   bool &write_input_data, int &f_precision, int &n_threads, double &threshold, int &times, char& delimiter);





void help();

void cmd_args(int argc, char* argv[]);
int cmd_args_v2(int argc, char** argv);
int cmd_args_v3(int argc, char** argv);

void version();

void version_notified(int a);

#endif //CONVOLUTION_CMD_ARGS_H
