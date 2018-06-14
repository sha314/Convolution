#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <map>

#include "src/binomial.h"
#include "src/convolution.h"
#include "src/data_reader.h"

using namespace std;


void time_complexity_data(){

    ofstream fout("Time complexity.txt");
    fout << "<N>\t<time version 1>\t<time version 2>\t<time version 3>" << endl;

    clock_t t_v1, t_v2;

    //generating random data
    for(size_t N=10000; N != 100000; N += 5000){
        cout << "N= " << N << endl;
        fout << N ;
        vector<double> data(N);
        for(size_t i{}; i != N; ++ i){
            data[i] = rand() / double(RAND_MAX);
        }
        t_v1 = clock();
        convolve_v1(data);
        fout << '\t' << (clock() - t_v1)/double(CLOCKS_PER_SEC) ;

        t_v2 = clock();
        convolve_v2(data);
        fout << '\t' << (clock() - t_v2)/double(CLOCKS_PER_SEC);

        t_v2 = clock();
        convolve_v3(data);
        fout << '\t' << (clock() - t_v2)/double(CLOCKS_PER_SEC);

        fout << endl;
    }

    fout.close();
}


void time_complexity_v3_data(){

    ofstream fout("Time complexity_multi.txt");
    fout << "#<N>\t<time version 1>\t<time version 2>" << endl;

    clock_t t_v1, t_v2;
    size_t number_of_columns = 3;
    //generating random data
    for(size_t N=10000; N != 20000; N += 500){
        cout << "N= " << N << endl;
        fout << N << '\t';
        vector<vector<double>> data(N);
        for(size_t i{}; i != N; ++ i){
            data[i] = vector<double>(number_of_columns);
            for(size_t j{}; j != data[i].size(); ++ j) {
                data[i][j] = rand() / double(RAND_MAX);
            }
        }
        t_v1 = clock();
        convolve_multi_v1(data);
        fout << (clock() - t_v1)/double(CLOCKS_PER_SEC);

        t_v2 = clock();
        convolve_multi_v2(data);
        fout << '\t' << (clock() - t_v2)/double(CLOCKS_PER_SEC);

        fout << endl;
    }

    fout.close();
}

template <typename T>
void print_vector(const std::vector<T>& vec){
    cout << '{';
    for(auto d: vec){
        cout << d << ',';
    }
    cout << '}' << endl;
}


void run_program(){
    cout << "Enter filename : ";
    string filename;
    getline(cin, filename);

    cout << "Enter independent column : ";
    unsigned ic{};
    string tmp;
    getline(cin, tmp);
    istringstream iss0(tmp);
    iss0 >> ic;
    cout << "Enter independent column name : ";
    string icolumn_name;
    getline(cin, icolumn_name);

    cout << "Enter columns to be convolved : ";
    string cols;
    getline(cin, cols);
    istringstream iss(cols);
    unsigned c;
    vector<unsigned> usecols;
    while(iss >> c){
        usecols.push_back(c);
    }

    cout << "Enter names of the columns to be convolved : ";
    string col_names, col_name;
    getline(cin, col_names);
    istringstream iss2(col_names);
    vector<string> usecols_names;
    while(iss2 >> col_name){
        usecols_names.push_back(col_name);
    }
    // display the given input

    cout << filename << endl;
    cout << ic << endl;
    cout << icolumn_name << endl;
    print_vector(usecols);
    print_vector(usecols_names);

    // run program from here
    // read convolution data from file
    map<string, unsigned> header = read_header(filename);
    vector<double> independent_data = loadtxt(filename, ic, header["data_line"]-1);
    vector<double> data = loadtxt(filename, usecols[0], header["data_line"]-1);
//    print_vector(data);
    cout << "len(data) " << data.size() << endl;
    // perform convolution
    vector<double> convolved_data = convolve_v1(data);

    // write independent column, original columns, convolved columns
    string outfilename = filename.substr(0, filename.find("_2018"));
    outfilename += "_convolved.txt";
    ofstream fout(outfilename);
    // with their names and appropriate header
    string new_header_info = "#";

    new_header_info += usecols_names[0];
    for(size_t i{1}; i != usecols_names.size(); ++i){
        new_header_info += "\t<" + usecols_names[i] + ">";
        new_header_info += "\t<" + usecols_names[i] + " convolved>";
    }

    fout << new_header_info << endl;
    fout << "#Convolved data" << endl;
    fout << "BEGIN_HEADER" << endl;
    fout << "ensemble_size\t" << header["ensemble_size"] << endl;
    fout << "length\t" << header["length"] << endl;
    fout << "data_line\t" << 7 << endl;
    fout << "END_HEADER" << endl;

    for(size_t i{}; i != independent_data.size(); ++i){
        fout << independent_data[i] << '\t' << data[i] << '\t' << convolved_data[i] << endl;
    }
    fout.close();
}



void run_program_multi(){
    cout << "Enter filename : ";
    string filename;
    getline(cin, filename);

    cout << "Enter independent column : ";
    unsigned ic{};
    string tmp;
    getline(cin, tmp);
    istringstream iss0(tmp);
    iss0 >> ic;
    cout << "Enter independent column name : ";
    string icolumn_name;
    getline(cin, icolumn_name);

    cout << "Enter columns to be convolved : ";
    string cols;
    getline(cin, cols);
    istringstream iss(cols);
    unsigned c;
    vector<unsigned> usecols;
    while(iss >> c){
        usecols.push_back(c);
    }

    cout << "Enter names of the columns to be convolved : ";
    string col_names, col_name;
    getline(cin, col_names);
    istringstream iss2(col_names);
    vector<string> usecols_names;
    while(iss2 >> col_name){
        usecols_names.push_back(col_name);
    }
    // display the given input

    cout << filename << endl;
    cout << ic << endl;
    cout << icolumn_name << endl;
    print_vector(usecols);
    print_vector(usecols_names);

    // run program from here
    // read convolution data from file
    cout << "Enter header type (0 for raw and 1 for JSON) : ";
    int header_type{-1};
    cin >> header_type;
    map<string, unsigned> header;
    vector<double> independent_data ;
    vector<vector<double>> data ;

    if(header_type == 0) {
        header = read_header(filename);
        independent_data = loadtxt(filename, ic, header["data_line"]-1);
        data = loadtxt(filename, usecols, header["data_line"]-1);
    }
    else if(header_type == 1){
        header = read_header_json(filename);
        cout << "Every line except header and data should be commented : line " << __LINE__ << endl;
        independent_data = loadtxt(filename, ic, 1);
        data = loadtxt(filename, usecols, 1);
    }else{
        cout << "Invalid header type" << endl;
    }
    cout << "Header info" << endl;
    for(auto i: header){
        cout << i.first << "=>" << i.second << endl;
    }

//    print_vector(independent_data);
//    for(auto v : data) {
//        print_vector(v);
//    }
//    cout << "len(data) " << data.size() << endl;
    // perform convolution
    vector<vector<double>> convolved_data = convolve_multi_v1(data);

    // write independent column, original columns, convolved columns
    string outfilename = filename.substr(0, filename.find("_2018"));
    outfilename += "_convolved.txt";
    ofstream fout(outfilename);
    // with their names and appropriate header

    if(header_type == 0) {
        fout << output_header_raw(icolumn_name, usecols_names, header) << endl;
    }
    else if(header_type == 1) {
        fout << output_header_json(icolumn_name, usecols_names, header) << endl;
    }
    for(size_t i{}; i != independent_data.size(); ++i){
        fout << independent_data[i];
        for(size_t j{}; j != data[i].size(); ++j){
            fout << '\t' << data[i][j] << '\t' << convolved_data[i][j];
        }
        fout << endl;
    }
    fout.close();
}


void test_run_program(){
    string filename = "sq_lattice_site_percolation_100_calculated_2018.1.31_21.21.3.txt" ;
    unsigned ic{0};
    vector<unsigned> usecols{4};


    // run program from here
    // read convolution data from file
    map<string, unsigned> header = read_header(filename);
    vector<double> independent_data = loadtxt(filename, ic, header["data_line"]-1);
    vector<double> data = loadtxt(filename, usecols[0], header["data_line"]-1);
//    print_vector(data);
    cout << "len(data) " << data.size() << endl;
    // perform convolution
    vector<double> convolved_data = convolve_v1(data);
    // write independent column, original columns, convolved columns
    string outfilename = filename.substr(0, filename.find("_2018"));
    outfilename += "_convolved.txt";
    ofstream fout(outfilename);
    // with their names and appropriate header
    string new_header_info = "#convolved data";
    fout << "#<p>\t<C>\t<colvolvec C>" << endl;
    for(size_t i{}; i != independent_data.size(); ++i){
        fout << independent_data[i] << '\t' << data[i] << '\t' << convolved_data[i] << endl;
    }
    fout.close();
}


void test_multi_run_program(){
    string filename = "sq_lattice_site_percolation_100_calculated_2018.1.31_21.21.3.txt" ;
    unsigned ic{0};
    vector<unsigned> usecols{4, 5};


    // run program from here
    // read convolution data from file
    map<string, unsigned> header = read_header(filename);
    vector<double> independent_data = loadtxt(filename, ic, header["data_line"]-1);
    vector<vector<double>> data = loadtxt(filename, usecols, header["data_line"]-1);
//    for(auto v : data) {
//        print_vector(v);
//    }
//    cout << "len(data) " << data.size() << endl;
    // perform convolution
    vector<vector<double>> convolved_data = convolve_multi_v1(data);
    // write independent column, original columns, convolved columns
    string outfilename = filename.substr(0, filename.find("_2018"));
    outfilename += "_convolved.txt";
    ofstream fout(outfilename);
    // with their names and appropriate header
    string new_header_info = "#";

//    new_header_info += usecols_names[0];
//    for(size_t i{0}; i != usecols.size(); ++i){
//        new_header_info += "\t<" + usecols_names[i+1] + ">";
//        new_header_info += "\t<" + usecols_names[i+1] + " convolved>";
//    }

    fout << new_header_info << endl;
    fout << "#Convolved data" << endl;
    fout << "BEGIN_HEADER" << endl;
    fout << "ensemble_size\t" << header["ensemble_size"] << endl;
    fout << "length\t" << header["length"] << endl;
    fout << "data_line\t" << 7 << endl;
    fout << "END_HEADER" << endl;
    for(size_t i{}; i != independent_data.size(); ++i){
        fout << independent_data[i];
        for(size_t j{}; j != data[i].size(); ++j){
            fout << '\t' << data[i][j] << '\t' << convolved_data[i][j];
        }
        fout << endl;
    }
    fout.close();
}


void cmd_args(int argc, char* argv[]){

}


void time_test(){
    time_complexity_data();
//    time_complexity_v3_data();

}


/***
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
    cout << "Convolution on big data" << endl;
    clock_t t = clock();

//    time_test();
//    test_run_program();
//    test_multi_run_program();

//    run_program();
    run_program_multi(); // this is it


//    print_vector(binomial_distribution_v1(10, 3, 0.5));

    cout << "Program finished " << endl;
    cout << "Time elapsed " << (clock() - t) / (double(CLOCKS_PER_SEC) * 60) << " minutes" << endl;
    return 0;
}

