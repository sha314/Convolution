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


#endif //CONVOLUTION_CONVOLUTION_H
