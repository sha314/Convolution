
#include "../convolution/convolution.h"
#include "../io/data_reader.h"
#include "test2.h"
#include <iostream>
#include <fstream>

using namespace std;


void test_factors(){
    Convolution convolution(100);
    auto f1 = convolution.factor_forward();
    auto f2 = convolution.factor_backward();
    for(size_t i{}; i < f1.size(); ++i){
        cout << i << " => " << f1[i] << ", " << f2[i] << endl;
    }
}

void test_convolution_basic(size_t N){

    Convolution convolution;

    vector<double> b_data_in(N);
    for(size_t i{}; i < N; ++i){
        b_data_in[i] = 1.0;
    }
    auto b_data_out = convolution.run(b_data_in);

    for(size_t i{}; i < b_data_in.size(); ++i){
        cout << b_data_in[i] << ", " << b_data_out[i] << endl;
    }

}

void test1_convolution(){
    string filename="data_json.txt";
    vector<double> a = loadtxt(filename, 0, 0, ' ');
    vector<double> b_data_in = loadtxt(filename, 1, 0, ' ');
    Convolution convolution;
    auto b_data_out = convolution.run(b_data_in);
    auto f_factor = convolution.factor_forward();
    auto b_factor = convolution.factor_backward();
    view(f_factor);
    view(b_factor);
//    ofstream fout(filename+"_convoluted.txt");
//    for(size_t i{}; i < b_data_in.size(); ++i){
//        fout << a[i] << '\t' << b_data_out[i] << endl;
//    }
//    fout.close();
}