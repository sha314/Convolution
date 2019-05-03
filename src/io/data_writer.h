//
// Created by shahnoor on 10/27/18.
//

#ifndef CONVOLUTION_DATA_WRITER_H
#define CONVOLUTION_DATA_WRITER_H

#include <iostream>
#include <string>
#include <vector>

void
savetxt_multi(
        const std::string &in_filename,
        const std::string &out_filename,
        const std::string &info,
        bool write_header_and_comment,
        char delimeter,
        bool write_input_data,
        const std::vector<std::vector<double>> &a_data,
        const std::vector<std::vector<double>> &b_data_in,
        const std::vector<std::vector<double>> &b_data_out,
        int precision
);



#endif //CONVOLUTION_DATA_WRITER_H
