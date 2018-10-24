//
// Created by shahnoor on 8/6/18.
//

#ifndef CONVOLUTION_PRINTER_H
#define CONVOLUTION_PRINTER_H

#include <vector>
#include <iostream>

template <typename T>
void print_vector(const std::vector<T>& vec){
    std::cout << '{';
    for(auto d: vec){
        std::cout << d << ',';
    }
    std::cout << '}' << std::endl;
}




#endif //CONVOLUTION_PRINTER_H
