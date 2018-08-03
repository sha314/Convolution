#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <map>
#include <chrono>

#include "src/binomial.h"
#include "src/convolution.h"
#include "src/data_reader.h"


void run_in_main();

using namespace std;


template <typename T>
void print_vector(const std::vector<T>& vec){
    cout << '{';
    for(auto d: vec){
        cout << d << ',';
    }
    cout << '}' << endl;
}



void take_input_cmd(string &filename, unsigned int &ic, string &icolumn_name, vector<unsigned int> &usecols,
                    vector<string> &usecols_names, int &header_type) {
    cout << "Enter filename : ";
    getline(cin, filename);

    cout << "Enter independent column : ";
    string tmp;
    getline(cin, tmp);
    istringstream iss0(tmp);
    iss0 >> ic;
    cout << "Enter independent column name : ";
    getline(cin, icolumn_name);

    cout << "Enter columns to be convolved : ";
    string cols;
    getline(cin, cols);
    istringstream iss(cols);
    unsigned c;
    while(iss >> c){
        usecols.push_back(c);
    }

    cout << "Enter names of the columns to be convolved : ";
    string col_names, col_name;
    getline(cin, col_names);
    istringstream iss2(col_names);
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
    cin >> header_type;
}



void take_input_json_cmd(string &filename, unsigned int &ic, string &icolumn_name, vector<unsigned int> &usecols,
                    vector<string> &usecols_names) {
    cout << "Enter filename : ";
    getline(cin, filename);

    cout << "Enter independent column : ";
    string tmp;
    getline(cin, tmp);
    istringstream iss0(tmp);
    iss0 >> ic;
    cout << "Enter independent column name : ";
    getline(cin, icolumn_name);

    cout << "Enter columns to be convolved : ";
    string cols;
    getline(cin, cols);
    istringstream iss(cols);
    unsigned c;
    while(iss >> c){
        usecols.push_back(c);
    }

    cout << "Enter names of the columns to be convolved : ";
    string col_names, col_name;
    getline(cin, col_names);
    istringstream iss2(col_names);
    while(iss2 >> col_name){
        usecols_names.push_back(col_name);
    }
    // display the given input

    cout << filename << endl;
    cout << ic << endl;
    cout << icolumn_name << endl;
    print_vector(usecols);
    print_vector(usecols_names);

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


/**
 *
 * @return Time (sec) required to calculate
 */
double run_program_multi(){
    string filename;
    unsigned int ic;
    string icolumn_name;
    vector<unsigned int> usecols;
    vector<string> usecols_names;
    int header_type;
    take_input_cmd(filename, ic, icolumn_name, usecols, usecols_names, header_type);

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
    auto t0 = std::chrono::system_clock::now();
    // perform convolution
//    vector<vector<double>> convolved_data = convolve_multi_v1(data); // single thread
//    vector<vector<double>> convolved_data = convolve_multi_threaded_v1(data); //multi thread
    vector<vector<double>> convolved_data = convolve_multi_threaded_v2(data); //multi thread

    auto t1 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t1-t0;

    cout << "calculation completed successfully : line " << __LINE__ << endl;

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

    return elapsed_seconds.count();
}


/**
 * Only JSON formated header is allowed
 * @return Time (sec) required to calculate
 */
double run_program_json_multi(string &filename, unsigned int &ic, string &icolumn_name, vector<unsigned int> &usecols,
                              vector<string> &usecols_names){

    map<string, unsigned> header;
    vector<double> independent_data ;
    vector<vector<double>> data ;

    header = read_header_json(filename);
    cout << "Every line except header and data should be commented : line " << __LINE__ << endl;
    independent_data = loadtxt(filename, ic, 1);
    data = loadtxt(filename, usecols, 1);

    cout << "Header info" << endl;
    for(auto i: header){
        cout << i.first << "=>" << i.second << endl;
    }


//    print_vector(independent_data);
//    for(auto v : data) {
//        print_vector(v);
//    }
//    cout << "len(data) " << data.size() << endl;
    auto t0 = std::chrono::system_clock::now();
    // perform convolution
//    vector<vector<double>> convolved_data = convolve_multi_v1(data); // single thread
//    vector<vector<double>> convolved_data = convolve_multi_threaded_v1(data); //multi thread
    vector<vector<double>> convolved_data = convolve_multi_threaded_v2(data); //multi thread

    auto t1 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t1-t0;

    cout << "calculation completed successfully : line " << __LINE__ << endl;

    // write independent column, original columns, convolved columns
    string outfilename = filename.substr(0, filename.find("_2018"));
    outfilename += "_convolved.txt";
    ofstream fout(outfilename);
    // with their names and appropriate header

    fout << output_header_json(icolumn_name, usecols_names, header) << endl;

    for(size_t i{}; i != independent_data.size(); ++i){
        fout << independent_data[i];
        for(size_t j{}; j != data[i].size(); ++j){
            fout << '\t' << data[i][j] << '\t' << convolved_data[i][j];
        }
        fout << endl;
    }
    fout.close();

    return elapsed_seconds.count();
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

void
parse_input_json_cmd(
        int argc,
        char* argv[],
        string &filename,
        unsigned int &ic,
        string &icolumn_name,
        vector<unsigned int> &usecols,
        vector<string> &usecols_names
){
    cout << "parse cmd inputs : line " << __LINE__ << endl;
    char *chs;
    for(int i{1}; i < argc;){
        chs = argv[i];
        switch(chs[1]){
            case 'f':
                ++i;
                filename = argv[i];
            case 'i':
                ++i;


        }
    }
}


void run_in_main(int argc, char* argv[]) {//    time_test();
//    test_run_program();
//    test_multi_run_program();

//    run_program();
//    double required_time = run_program_multi(); // this is it

    string filename;
    unsigned int ic;
    string icolumn_name;
    vector<unsigned int> usecols;
    vector<string> usecols_names;

    take_input_json_cmd(filename, ic, icolumn_name, usecols, usecols_names);
    double required_time = run_program_json_multi(filename, ic, icolumn_name, usecols, usecols_names);


//    print_vector(binomial_distribution_v1(10, 3, 0.5));

    cout << "Time elapsed for convolution only " << required_time / 60.0 << " minutes" << endl;
}


void test_in_main(int argc, char **argv){
    string filename = "sq_lattice_site_percolation_periodic_200-calculated-.txt";
    vector<double> independent_data = loadtxt(filename, 0, 1);

    vector<double> data = loadtxt(filename, 1, 1);

    cout << "data.size() " << data.size() << endl;


    Convolution conv;
    vector<double> out_data = conv.run_omp(data);
    conv.timeElapsed();

    string outfilename = filename + "_convolved.txt";
    ofstream fout(outfilename);
    // with their names and appropriate header


    for(size_t i{}; i != independent_data.size(); ++i){
        fout << independent_data[i] ;
        fout << '\t' << data[i];
        fout << '\t' << out_data[i];

        fout << endl;
    }
    fout.close();
}

void test_multi_in_main(int argc, char **argv){
    string filename = "sq_lattice_site_percolation_periodic_200-calculated-.txt";
    vector<double> independent_data = loadtxt(filename, 0, 1);
    vector<unsigned> cols = {1,2};
    vector<vector<double>> data = loadtxt(filename, cols, 1);

    cout << "data.size() " << data.size() << endl;
    cout << "data[0].size() " << data[0].size() << endl;

    for(size_t i{}; i != data[0].size(); ++i){
        for(size_t j{}; j != data.size(); ++j){
            cout << data[j][i] << ',';
            if(j > 5){
                break;
            }
        }
        cout << endl;
        if(i > 5){
            break;
        }
    }
    Convolution conv;
    vector<vector<double>> out_data = conv.run_multi_omp(data);
    conv.timeElapsed();

    string outfilename = filename + "_convolved.txt";
    ofstream fout(outfilename);
    // with their names and appropriate header


    for(size_t i{}; i != independent_data.size(); ++i){
        fout << independent_data[i] ;
        for(size_t j{}; j != data[0].size(); ++j) {
            fout << '\t' << data[i][j] << '\t' << out_data[i][j];
        }
        fout << endl;
    }
    fout.close();
}

/***
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
    cout << "Convolution of big data" << endl;

    auto t0 = std::chrono::system_clock::now();

//    run_in_main(argc, argv);
//    test_in_main(argc, argv);
    test_multi_in_main(argc, argv);

    auto t1 = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = t1-t0;
    std::time_t end_time = std::chrono::system_clock::to_time_t(t1);


    cout << "Program finished at " << std::ctime(&end_time) << endl;
    cout << "Total Time elapsed " << elapsed_seconds.count()/60 << " minutes" << endl;
    return 0;
}

