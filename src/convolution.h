//
// Created by shahnoor on 2/1/2018.
//

#ifndef CONVOLUTION_CONVOLUTION_H
#define CONVOLUTION_CONVOLUTION_H

#include <vector>
#include <cstddef>

std::vector<double> convolve_v1(const std::vector<double>& data);
std::vector<double> convolve_v2(const std::vector<double>& data);
std::vector<double> convolve_v3(const std::vector<double>& data);

std::vector<std::vector<double>> convolve_multi_v1(const std::vector<std::vector<double>>& data);
std::vector<std::vector<double>> convolve_multi_v2(const std::vector<std::vector<double>>& data);
// 7 July, 2018
std::vector<std::vector<double>> convolve_multi_threaded_v1(const std::vector<std::vector<double>> &data_in);
// 8 July, 2018
std::vector<std::vector<double>> convolve_multi_threaded_v2(const std::vector<std::vector<double>> &data_in);

// 7 July, 2018
void calculate_convolution_for_rows(
        const size_t start,
        const size_t end,
        const std::vector<std::vector<double>> &data,
        std::vector<std::vector<double>> &out_data
);

// 8 july, 2018
void calculate_convolution_for_rows_v2(
        const size_t start,
        const size_t end,
        const std::vector<std::vector<double>> &data,
        std::vector<std::vector<double>> &out_data,
        const size_t id,
        std::vector<size_t> &iterations
);


#endif //CONVOLUTION_CONVOLUTION_H
