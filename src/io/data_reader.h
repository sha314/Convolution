//
// Created by shahnoor on 2/2/2018.
//

#ifndef CONVOLUTION_DATA_READER_H
#define CONVOLUTION_DATA_READER_H


#include <string>
#include <vector>
#include <map>

std::map<std::string, unsigned> read_header(std::string filename, char delemiter=' ', char comment='#');
std::map<std::string, unsigned> read_header_json(std::string filename, char comment='#');

std::vector<double> loadtxt(std::string filename, int usecols,
                       int skiprows, char delemiter=' ', char comment='#');

std::vector<std::vector<double>> loadtxt(std::string filename, const std::vector<int>& usecols,
                                         int skiprows, char delemiter=' ', char comment='#');

std::vector<std::vector<double>> loadtxt_v2(std::string filename, const std::vector<int>& usecols,
                                         int skiprows, char delemiter=' ', char comment='#');


std::vector<std::string> explode_to_string(const std::string &str, const char &ch);
std::vector<int>         explode_to_int(const std::string &str, const char &ch);
std::vector<double>      explode_to_float(const std::string &s, const char &c);

std::string output_header_json(
        const std::string& icolumn_name,
        const std::vector<std::string>& usecols_names,
        const std::map<std::string, unsigned>& header);

std::string output_header_raw(
        const std::string& icolumn_name,
        const std::vector<std::string>& usecols_names,
        std::map<std::string, unsigned>& header
);

#endif //CONVOLUTION_DATA_READER_H
