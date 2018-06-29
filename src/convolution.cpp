//
// Created by shahnoor on 2/1/2018.
//

#include <cmath>
#include <iostream>
#include "convolution.h"
#include "binomial.h"

using namespace std;


/**
 * Perform convolution of a given data list (single column) using the following formula
 *    X_p = B(N, n, p) * X_n
 *    where, B(N, n, p) = N^C_n * p^n * (1-p)^(N-n)
 * Since it is costly to use this formula directly, this function uses a simplified one
 *    B(N,i,p) = 1  # is set arbitrarily
 *    B(N,n,p) = B(N, n-1, p) * ((N-n+1)/n) * (p/(1-p))    # for i > n
 *    B(N,n,p) = B(N, n+1, p) * (n+1)/(N-n)) * ((1-p)/p)   # for i < n
 *
 * Complexity, O(2.426)
 *
 * @param data : array of double valued data
 * @return     : array of double valued convolved data
 */
vector<double> convolve_v1(const vector<double>& data) {

    size_t size = data.size();
    vector<double> binomial;
    vector<double> out_data(size);
    double p;
//    size_t len = size / 100;
    for(size_t i{}; i != size; ++i) {
//        if(i % len == 0){
//            cout << "Running for row " << i << endl;
//        }
        p = (i + 1) / double(size);
        binomial = binomial_distribution_v2(size, i, p);
        for(size_t j{}; j != size; ++j){
            out_data[i] += binomial[j]*data[j];
        }
    }
    return out_data;
}



/**
 * Perform convolution of a given data list (single column) using the following formula
 *    X_p = B(N, n, p) * X_n
 *    where, B(N, n, p) = N^C_n * p^n * (1-p)^(N-n)
 * Since it is costly to use this formula directly, this function uses a simplified one
 *    B(N,i,p) = 1  # is set arbitrarily
 *    B(N,n,p) = B(N, n-1, p) * ((N-n+1)/n) * (p/(1-p))    # for i > n
 *    B(N,n,p) = B(N, n+1, p) * (n+1)/(N-n)) * ((1-p)/p)   # for i < n
 *
 * Complexity, O(2.769)
 *
 * @param data : array of double valued data
 * @return     : array of double valued convolved data
 */
vector<double> convolve_v2(const vector<double>& data) {

    size_t size = data.size();
    vector<double> out_data(size);
    double probability, a, b, sm, normalization_const, factor;
//    size_t len = size / 100;
    for(size_t i{}; i != size; ++i) {
//        if(i % len == 0){
//            cout << "Running for row " << i << endl;
//        }

        probability = (i + 1) / double(size);
        a = b = 1;
        sm = data[i];   //# sum
        normalization_const = 1;
        for(size_t j{i+1}; j != size; ++j) {  // for j > i
            factor = (size - j + 1) * probability / (j * (1 - probability));
            b = a * factor;
            normalization_const += b;
            sm += b * data[j];
            a = b;  //# swapping
        }
        a = b = 1; //# reassign to unity
        for(size_t j{i-1}; j != -1; --j) {  // for j < i
            factor = (j + 1) * (1 - probability) / ((size - j) * probability);
            a = b * factor;
            sm += a * data[j];
            normalization_const += a;
            b = a;  //# swapping
        }
        sm /= normalization_const;  //# normalize    the result
        out_data[i] = sm;
    }
    return out_data;
}

/**
 * Perform convolution of a given data list (single column) using the following formula
 *    X_p = B(N, n, p) * X_n
 *    where, B(N, n, p) = N^C_n * p^n * (1-p)^(N-n)
 * Since it is costly to use this formula directly, this function uses a simplified one
 *    B(N,i,p) = 1  # is set arbitrarily
 *    B(N,n,p) = B(N, n-1, p) * ((N-n+1)/n) * (p/(1-p))    # for i > n
 *    B(N,n,p) = B(N, n+1, p) * (n+1)/(N-n)) * ((1-p)/p)   # for i < n
 *
 * Complexity, O(2.628)
 * todo Uses symmetry of the binomial ??
 * @param data : array of double valued data
 * @return     : array of double valued convolved data
 */
vector<double> convolve_v3(const vector<double>& data) {

    size_t size = data.size();
    vector<double> binomial(size);
    vector<double> out_data(size);
    double p;
//    size_t len = size / 100;
    for(size_t i{}; i != size; ++i) {
//        if(i % len == 0){
//            cout << "Running for row " << i << endl;
//        }
        p = (i + 1) / double(size);
        binomial[i] = 1;
        double factor;
        double normalization_const={1};
        for (size_t j{i + 1}; j != size; ++j) {  // for j > center
            factor = (size - j + 1) * p / (j * (1 - p));
            binomial[j] = binomial[j - 1] * factor;
            normalization_const += binomial[j];
        }

        for (size_t j{i - 1}; j != -1; --j) {  // for j < center
            factor = (j + 1) * (1 - p) / ((size - j) * p);
            binomial[j] = binomial[j + 1] * factor;
            normalization_const += binomial[j];
        }

        for(size_t j{}; j != size; ++j){
            out_data[i] += binomial[j]*data[j]/normalization_const;
        }
    }
    return out_data;
}


/**
 * Perform convolution of a given data list (multiple column) using the following formula
 *    X_p = B(N, n, p) * X_n
 *    where, B(N, n, p) = N^C_n * p^n * (1-p)^(N-n)
 * Since it is costly to use this formula directly, this function uses a simplified one
 *    B(N,i,p) = 1  # is set arbitrarily
 *    B(N,n,p) = B(N, n-1, p) * ((N-n+1)/n) * (p/(1-p))    # for i > n
 *    B(N,n,p) = B(N, n+1, p) * (n+1)/(N-n)) * ((1-p)/p)   # for i < n
 *
 * Complexity, O(2.482) with 3 columns
 *
 * @param data : n-dimensional array of double valued data
 *      for example: following data has 3 columns and N rows
 *              0.25    0.454   0.548
 *              0.457   0.187   0.154
 *              0.578   0.951   0.487
 *              .       .       .
 *              .       .       .
 *              .       .       .
 * @return     : n-dimensional array of double valued convolved data
 */
vector<vector<double>> convolve_multi_v1(const vector<vector<double>>& data) {
    size_t size = data.size();
    vector<vector<double>> out_data(data.size());
    vector<double> binomial;
    double probability;
    size_t len = size / 1000;

    for(size_t i{}; i != size; ++i) {
        if(i % len == 0){
//          cout << "Running for row " << i << endl;
//            cout << "\33[2K"; // erase the current line
            cout << '\r'; // return the cursor to the start of the line
            cout << "at " << 100 * i / double(size) << "%";
            std::fflush(stdout);
//            cout << endl;
        }

        probability = (i + 1) / double(size);
        binomial = binomial_distribution_v2(size, i, probability);

        out_data[i] = vector<double>(data[0].size());
        for (size_t k{}; k != data[0].size(); ++k) {
            for (size_t j{}; j != size; ++j) {
                out_data[i][k] += data[j][k] * binomial[j];
            }
        }
    }
    return out_data;
}



/**
 * Perform convolution of a given data list (multiple column) using the following formula
 *    X_p = B(N, n, p) * X_n
 *    where, B(N, n, p) = N^C_n * p^n * (1-p)^(N-n)
 * Since it is costly to use this formula directly, this function uses a simplified one
 *    B(N,i,p) = 1  # is set arbitrarily
 *    B(N,n,p) = B(N, n-1, p) * ((N-n+1)/n) * (p/(1-p))    # for i > n
 *    B(N,n,p) = B(N, n+1, p) * (n+1)/(N-n)) * ((1-p)/p)   # for i < n
 *
 * Complexity, O(3.184) with 3 columns
 *
 * @param data : n-dimensional array of double valued data
 *      for example: following data has 3 columns and N rows
 *              0.25    0.454   0.548
 *              0.457   0.187   0.154
 *              0.578   0.951   0.487
 *              .       .       .
 *              .       .       .
 *              .       .       .
 * @return     : n-dimensional array of double valued convolved data
 */
vector<vector<double>> convolve_multi_v2(const vector<vector<double>>& data) {
    size_t size = data.size();
    vector<vector<double>> out_data(data.size());
    vector<double> sum(data[0].size());
    vector<double> binomial;
    double probability, a, b, normalization_const, factor;
//    size_t len = size / 100;
    for(size_t i{}; i != size; ++i) {
//        if(i % len == 0){
//            cout << "Running for row " << i << endl;
//        }

        probability = (i + 1) / double(size);
        a = b = 1;
        for(size_t k{}; k != sum.size(); ++k) {
            sum[k] = data[i][k];   //# sum
        }
        normalization_const = 1;
        for(size_t j{i+1}; j != size; ++j) {  // for j > i
            factor = (size - j + 1) * probability / (j * (1 - probability));
            b = a * factor;
            normalization_const += b;
            for(size_t k{}; k != sum.size(); ++k) {
                sum[k] += b*data[j][k];   //# sum
            }
            a = b;  //# swapping
        }
        a = b = 1; //# reassign to unity
        for(size_t j{i-1}; j != -1; --j) {  // for j < i
            factor = (j + 1) * (1 - probability) / ((size - j) * probability);
            a = b * factor;
            for(size_t k{}; k != sum.size(); ++k) {
                sum[k] += a*data[j][k];   //# sum
            }
            normalization_const += a;
            b = a;  //# swapping
        }
        for(size_t k{}; k != sum.size(); ++k) {
            sum[k] /= normalization_const;   //normalize the result
        }

        out_data[i] = sum;
    }
    return out_data;
}

