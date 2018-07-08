//
// Created by shahnoor on 2/1/2018.
//

#include <cmath>
#include <iostream>
#include <thread>
#include <mutex>
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
 * Progress information in the console.
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
            cout << "\33[2K"; // erase the current line
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



/**
 *
 * @param data
 * @param out_data
 * @param i
 */
void calculate_convolution_for_rows(
        const size_t start,
        const size_t end,
        const vector<vector<double>> &data,
        vector<vector<double>> &out_data
) {
    size_t size = data.size();
    double probability{};
    vector<double> binomial;
    double sum{};
    for(size_t i{start}; i != end; ++i) {
        probability = (i + 1) / double(size);
        binomial = binomial_distribution_v2(size, i, probability);
        for (size_t k{}; k != data[0].size(); ++k) {
            sum = 0;
            for (size_t j{0}; j != size; ++j) {
                 sum += data[j][k] * binomial[j];
            }
            out_data[i][k] = sum;
        }
    }
}



/**
 *
 * @param data
 * @param out_data
 * @param i
 */
void calculate_convolution_for_rows_v2(
        const size_t start,
        const size_t end,
        const vector<vector<double>> &data,
        vector<vector<double>> &out_data,
        const size_t id,
        vector<size_t> &iterations
) {
    size_t size = data.size();
    double probability{};
    vector<double> binomial;
    double sum{};
    for(size_t i{start}; i != end; ++i) {
        // keeps track of the number of iteration  performed
        // so that it can be accessed from another thread
        ++iterations[id];

        probability = (i + 1) / double(size);
        binomial = binomial_distribution_v2(size, i, probability);
        for (size_t k{}; k != data[0].size(); ++k) {
            sum = 0;
            for (size_t j{0}; j != size; ++j) {
                sum += data[j][k] * binomial[j];
            }
            out_data[i][k] = sum;
        }
    }
}

/**
 * Multithreading is used. // todo
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
 * @param data_in : n-dimensional array of double valued data
 *      for example: following data has 3 columns and N rows
 *              0.25    0.454   0.548
 *              0.457   0.187   0.154
 *              0.578   0.951   0.487
 *              .       .       .
 *              .       .       .
 *              .       .       .
 * @return     : n-dimensional array of double valued convolved data
 */
vector<vector<double>> convolve_multi_threaded_v1(const vector<vector<double>> &data_in) {
    size_t size = data_in.size();
    vector<vector<double>> data_out(data_in.size());

    for(size_t i{}; i != size; ++i) { // initialize sizes
        data_out[i] = vector<double>(data_in[0].size());
    }

    size_t len = size / 1000;

    size_t number_of_threads = std::thread::hardware_concurrency(); // number of threads

    size_t division = size / number_of_threads;

    vector<thread> threads(number_of_threads);

    size_t start, end;
    for(size_t i{}; i != number_of_threads; ++i){
        start = i*division;
        end = (i+1) * division;
        // each thread finishes independently
        threads[i] = std::thread(calculate_convolution_for_rows, start, end, std::ref(data_in), std::ref(data_out));

    }

    // join the threads
    for(size_t j{}; j != number_of_threads; ++j){
        if(threads[j].joinable())  threads[j].join();
    }

    return data_out;
}


/**
 * Multithreading is used with progress information in the console.
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
 * @param data_in : n-dimensional array of double valued data
 *      for example: following data has 3 columns and N rows
 *              0.25    0.454   0.548
 *              0.457   0.187   0.154
 *              0.578   0.951   0.487
 *              .       .       .
 *              .       .       .
 *              .       .       .
 * @return     : n-dimensional array of double valued convolved data
 */
vector<vector<double>> convolve_multi_threaded_v2(const vector<vector<double>> &data_in) {
    size_t size = data_in.size();
    vector<vector<double>> data_out(data_in.size());

    for(size_t i{}; i != size; ++i) { // initialize sizes
        data_out[i] = vector<double>(data_in[0].size());
    }

    size_t number_of_threads = std::thread::hardware_concurrency(); // number of threads

    size_t division = size / number_of_threads;

    vector<thread> threads(number_of_threads); // holds the threads
    vector<size_t> iterations(number_of_threads); // progress report

    size_t start, end;

    std::mutex mutex_main_thread;
    for(size_t i{}; i != number_of_threads; ++i){
        start = i*division;
        end = (i+1) * division;
        // each thread finishes independently
        threads[i] = std::thread(
                calculate_convolution_for_rows_v2,
                start,
                end,
                std::ref(data_in),
                std::ref(data_out),
                i,
                std::ref(iterations)
        );

    }


    // summary in the main thread
    size_t sum{};
    while(sum < size)
    {
        sum = 0;
        {
            std::lock_guard<mutex> _locker(mutex_main_thread);
            for (size_t i{}; i != iterations.size(); ++i) {
                sum += iterations[i];
            }
        }
        cout << "\33[2K"; // erase the current line
        cout << '\r'; // return the cursor to the start of the line
        cout << "progress : " << 100 * sum / double(size)  << "%";
        std::fflush(stdout);
        std::this_thread::sleep_for(std::chrono::duration<double>(5));
    }

    cout << endl;
    // join the threads
    for(size_t j{}; j != number_of_threads; ++j){
        if(threads[j].joinable())  threads[j].join();
    }

    return data_out;
}