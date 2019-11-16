//
// Created by shahnoor on 4/14/19.
//

#ifndef CONVOLUTION_TEST2_H
#define CONVOLUTION_TEST2_H

#include<iostream>
#include<vector>
#include<cmath>
#include <typeinfo>

using namespace std;


std::vector<double> factorials(size_t N);

std::vector<long double> nCr_list(size_t N);

double nCr(size_t N, size_t r, const std::vector<double>& facts);

template <class T>
void view(std::vector<T> & arr){
    std::cout << "viewing : line " << __LINE__ << std::endl;
    std::cout << '{';
    for(auto a: arr){
        std::cout << a << ',';
    }
    std::cout << '}' << std::endl;
}

template <class T>
void view_matrix(std::vector<std::vector<T>>& arr){
    std::cout << "viewing : line " << __LINE__ << std::endl;
    std::cout << '{';
    for(auto a: arr){
        std::cout << '{';
        for(auto b: a) {
            std::cout << b << ',';
        }
        std::cout << '}' << std::endl;
    }
    std::cout << '}' << std::endl;
}

/**

	Data_size	Time
	100,000		11m46s
**/
std::vector<double> convolution(std::vector<double>& data_in);
std::vector<double> convolution_v2(std::vector<double>& data_in);


/**
Perform convolution and then derivative at once

	TODO
**/
std::vector<double> convolution_derivative(std::vector<double>& data_in);


/***
 * Test function
 *
 *
 */
void test2_convolution();
void test3_convolution();
void test4_convolution();

#endif //CONVOLUTION_TEST2_H
