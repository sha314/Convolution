//
// Created by shahnoor on 10/27/18.
//

#include "data_writer.h"
#include "../include/string_methods.h"
#include <iomanip>
#include <fstream>


using namespace std;

void
savetxt_multi(
        const string &in_filename,
        const string &out_filename,
        const string &info,
        bool write_header_and_comment,
        char delimeter,
        bool write_input_data,
        const vector<vector<double>> &a_data,
        const vector<vector<double>> &b_data_in,
        const vector<vector<double>> &b_data_out,
        int precision
) {
    ofstream fout(out_filename);
    if(write_header_and_comment) {
        ifstream fin(in_filename);
        string str;
        while (getline(fin, str)) {
            auto trimed = trim(str, ' ');
            if(isdigit(trimed[0])){
                break;
            }
            cout << str << endl;
            fout << str << endl;
        }
        fin.close();
    }
    fout << '#' << info << endl; // info cannot contain a new line character
    fout << "#convolved data" << endl;
    cout << b_data_out.size() << ", " << b_data_out[0].size() << endl;

    for(size_t i{}; i < b_data_in.size(); ++i){
        for(size_t j{}; j < a_data[0].size(); ++j){
            fout << setprecision(precision) << a_data[i][j];
        }
        for(size_t j{}; j < b_data_in[0].size(); ++j){
            if(write_input_data){
                fout << delimeter << setprecision(precision) << b_data_in[i][j];
            }
            fout << delimeter << setprecision(precision) << b_data_out[i][j];
//            cout << delimeter << setprecision(precision) << b_data_out[i][j];
        }
        fout << endl;
    }
    fout.close();
}