//
// Created by shahnoor on 4/14/19.
//

#include <fstream>
#include "test2.h"
#include "../io/data_reader.h"
#include "../convolution/convolution.h"
#include "../io/logger.h"

std::vector<double> factorials(size_t N){
    std::vector<double> fact(N+1);
    fact[0] = 1;
    for(size_t i{1}; i <= N; ++i){
        fact[i] = i * fact[i-1];
    }
    return fact;
}

std::vector<long double> nCr_list(size_t N){
    std::vector<long double> combinations(N+1);

    combinations[0] = 1;
    for(size_t i{1}; i <= N; ++i){
        combinations[i] = combinations[i-1] * (N-i+1)/i;
    }
    return combinations;
}



double nCr(size_t N, size_t r, const std::vector<double>& facts){
    return facts[N] / (facts[r] * facts[N-r]);
}


/**

	Data_size	Time
	100,000		11m46s
**/
std::vector<double> convolution(std::vector<double>& data_in){
    size_t N = data_in.size();
    auto binomial =  nCr_list(N+1); // all the combinations are here
    cout << "binomial.size() " << binomial.size() << endl;
    std::vector<double> data_out(N);
    double sum=0;

    for(size_t i{0}; i < N; ++i){
        double p= (i+1)/double(N);
        double one_mimus_p = 1-p;

        double f_p = 1;
        double prob= p/(1-p);
        double b_p= pow(one_mimus_p, N);

        sum = 0;
        for(size_t j{}; j <N; ++j){
            f_p *= prob;
            sum += binomial[j+1] * f_p * b_p * data_in[j];

        }
        data_out[i] = sum;
    }
    return data_out;
}

std::vector<double> convolution_v2(std::vector<double>& data_in){
    size_t N = data_in.size();
    auto binomial =  nCr_list(N+1); // all the combinations are here
//    view(binomial);
    cout << "binomial.size() " << binomial.size() << endl;
    std::vector<double> data_out(N);
    double sum=0;

    for(size_t i{0}; i < N; ++i){
        double p= (i+1)/double(N);

        double f_p = 1;
        double prob= p/(1-p);
        double constant = pow(1-p, N);
//        cout << "const " << constant << endl;
        sum = 0;
        for(size_t j{}; j <N; ++j){
            f_p *= prob;
//            if (f_p > 0) {
//                cout << "f_p " << f_p << endl;
//            }
//            cout << "binomial[j+1] " << binomial[j+1] << endl;
            sum += binomial[j+1] * f_p * data_in[j];

        }
        data_out[i] = sum * constant;
    }
    return data_out;
}

std::vector<double> convolution_v3(std::vector<double>& data_in){
    size_t N = data_in.size();
    vector<double> binom_factor(N);
    for(size_t i{1}; i <= N; ++i){
        // i cannot be zero
        binom_factor[i-1] = (N-i+1)/i;
    }
//    view(binomial);
    std::vector<double> data_out(N);
    double sum=0;

    for(size_t i{0}; i < N; ++i){
        double p= (i+1)/double(N);

        double f_p = 1;
        double prob= p/(1-p);
        double constant = pow(1-p, N);
//        cout << "const " << constant << endl;
        sum = 0;
        double factor = 1;
        for(size_t j{}; j <N; ++j){
            f_p *= prob;
//            if (f_p > 0) {
//                cout << "f_p " << f_p << endl;
//            }
//            cout << "binomial[j+1] " << binomial[j+1] << endl;
            factor *= binom_factor[j];
            sum += binom_factor[j] * f_p * data_in[j];

        }
        data_out[i] = sum * constant;
    }
    return data_out;
}




/**
Perform convolution and then derivative at once

	TODO
**/
std::vector<double> convolution_derivative(std::vector<double>& data_in){
    size_t N = data_in.size();
    auto binomial =  nCr_list(N+1); // all the combinations are here
    cout << binomial.size() << endl;
    std::vector<double> data_out(N);
    double sum=0;

    for(size_t i{0}; i < N; ++i){
        double p= (i+1)/double(N);

        double constant = pow(1-p, N-1)/p;
        double prob_factor = p/(1-p);
        double f_p = 1;

        sum = 0;
        for(size_t j{}; j <N; ++j){
            f_p *= prob_factor;
            sum += binomial[j+1] * f_p  * (j - N*p) * data_in[j];

        }
        data_out[i] = sum * constant;
    }
    return data_out;
}


void test2_convolution(){
    string filename="data_json.txt";
    vector<double> a = loadtxt(filename, 0, 3, ' ');
    vector<double> b_data_in = loadtxt(filename, 1, 3, ' ');
    cout << a[10] << ", " << b_data_in[10] << endl;
    auto b_data_out = convolution_v2(b_data_in);
    ofstream fout(filename+"_convoluted.txt");
    for(size_t i{}; i < b_data_in.size(); ++i){
        fout << a[i] << '\t' << b_data_out[i] << endl;
    }
    fout.close();
}

void test3_convolution(){
    string filename="data_json.txt";
    vector<double> a(10000, 1);

    Logging* log = Logging::getInstance();

    log->setFilename("e-9.txt");
    auto b_data_out = convolve_1d_fast(a, 1, 1e-9);

    log->setFilename("e-15.txt");
    convolve_1d_fast(a, 1, 1e-15);

    log->setFilename("e-18.txt");
    convolve_1d_fast(a, 1, 1e-18);

    log->setFilename("0.txt");
    convolve_1d_fast(a, 1, 0);

//    auto b_data_out = convolve_1d(a, 1);
//    ofstream fout(filename+"_convoluted.txt");
//    for(size_t i{}; i < b_data_out.size(); ++i){
//        fout << a[i] << '\t' << b_data_out[i] << endl;
//    }
//    fout.close();
}

void test4_convolution(){
    string filename="sq_lattice_site_percolation_periodic__400_2018.6.30_7.51.51.txt";
    filename ="sq_lattice_site_percolation_periodic_free-energy_L400_2019-11-13_205909.txt_convoluted.txt_derivative.txt";
    filename="sq_lattice_site_percolation_periodic_free-energy_L400_2019-11-13_205909.txt_convoluted.txt_derivative.txt_convoluted.txt";
    auto delimeter = analyze_delimeter(filename, 0, ' ');
    auto a = loadtxt_v2(filename,{1}, 0, delimeter);
    cout << a.size() << ", " << a[0].size() << endl;
    cout << "line " << __LINE__ << endl;
//    exit(0);
    auto b_data_out = convolve_2d_fast_diff(a, 1, 1, 1e-15);
//    auto b_data_out = convolve_1d(a, 1);
    cout << b_data_out.size() << ", " << b_data_out[0].size() << endl;
    ofstream fout(filename+"_convoluted.txt");
    for(size_t i{}; i < b_data_out.size(); ++i){
        for(size_t j{}; j < b_data_out[0].size(); ++j) {
            fout << a[i][j] << '\t' << b_data_out[i][j] << '\t';
//            cout << a[i][j] << '\t' << b_data_out[i][j] << '\t';
        }
        fout << endl;
//        cout << endl;
    }
    fout.close();
}
