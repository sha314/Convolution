//
// Created by shahnoor on 11/18/19.
//

#include "cmd_args.h"
#include "io/data_reader.h"
#include <iostream>
#include "boost/program_options.hpp"
#include "include/string_methods.h"
#include "convolution/convolution.h"
#include "io/data_writer.h"
#include "include/printer.h"

using namespace std;

void version(){
    cout << "version " << endl;
}

void version_notified(int a){
    cout << "version " << endl;
}


void get_option_b(int argc, char *const *argv, vector<int> &b_usecols, vector<string> &b_names, int i) {
    if(i < argc) {
        string tmp = argv[i];
        cout << tmp << endl;
        int sep = tmp.find(':');
        if(sep < 0){
            // no name is provided
            b_usecols = explode_to_int(tmp, ',');
            for(size_t k{}; k < b_usecols.size(); ++k){
                b_names.push_back("<>");
            }
        }else {
            string first = tmp.substr(0, sep);

            string second = tmp.substr(sep + 1);
//                        cout << "first " << first << endl;
//                        cout << "second " << second << endl;

            b_usecols = explode_to_int(first, ',');
            b_names = explode_to_string(second, ',');
            if(b_names.size() != b_usecols.size()){
                cout << "not enough column index or name" ;
                cout << ": line " << __LINE__ ; // comment this on deployment
                cout << endl;
            }
        }
    }
}

void get_option_a(int argc, char *const *argv, vector<int> &a_usecols, vector<string> &a_names, int i) {
    if(i < argc) {
        string tmp = argv[i];
        cout << tmp << endl;
        int sep = tmp.find(':');
        if(sep < 0){
            // no name is provided
            a_usecols = explode_to_int(tmp, ',');
            for(size_t k{}; k < a_usecols.size(); ++k){
                a_names.push_back("<>");
            }
        }else {
            string first = tmp.substr(0, sep);

            string second = tmp.substr(sep + 1);
//                        cout << "first " << first << endl;
//                        cout << "second " << second << endl;

            a_usecols = explode_to_int(first, ',');
            a_names = explode_to_string(second, ',');
            if(a_names.size() != a_usecols.size()){
                cout << "not enough column index or name" ;
                cout << ": line " << __LINE__ ; // comment this on deployment
                cout << endl;
            }
        }
    }
}


void help(){
    string hlp = R"***(Usage:
convolution [-f <STRING>] [-a <INT>,<INT>,...[:[<STRING>,<STRING>,...]]] [-b <INT>,<INT>,...[:[<STRING>,<STRING>,...]]] [-h] [-t <INT>] [-i <STRING>]

perform convolution based on provided options.

Options                      Description
  -a,                        columns that we want in the output file without performing convolution.
                             No default value.
  -b,                        columns that we want in the output file with performing convolution.
                             No default value.
  -c,                        If provided the header and comment from the input file will be written
                             without modification to the output file. Header is the first line of the
                             input file.
  -d                         Delimiter to use. Default value is ' '.
  -f                         name of the input file that we want to convolute. No default value.
  -i                         Info to write as comment in the output file
  -o                         name of the output file. If not provided the string '_convoluted.txt' will be
                             appended to the input file.
  -p, --precision            Floating point precision when writing in the data file. Default value is 10
  -s                         Number of rows to skip from the input file. Default value is 0.
  -t                         to test the performance of the convolution program. No default value.
      --threads              Explicitly specify number of thread to use. Default is the max number of thread
                             allowed by the system.
  -h, --help                 display this help and exit
  -v, --version              output version information and exit
  -w                         If provided input b data will be written to the output file.


The INT argument is an integer.
The STRING argument is a string of characters.

A line that begins with '#' is considered a commented line.

Exit status:
 0  if OK,
 1  if minor problems (e.g., cannot access subdirectory),
 2  if serious trouble (e.g., cannot access command-line argument).
)***";
    cout << hlp << endl;
}



void cmd_args(int argc, char* argv[]){
    cout << "Arguments" << endl;
    string in_filename, out_filename;
    string out_file_flag = "_convoluted.txt";
    int header_line{0};
    size_t test_size{0};
    vector<int> a_usecols, b_usecols; // b_usecols will be convolved and a_usecols will remain unchanged
    vector<string> a_names, b_names;
    string info;
    bool write_header_and_comment{false};
    int skiprows{0};
    char delimeter{' '};
    bool write_input_data {false};
    int flg;
    int f_precision{10};
    int n_threads{-1};
    if(argc == 1){
        help();
        exit(0);
    }
    for(int i{1}; i < argc;){
        cout << argv[i] << endl;

        flg = str2int(argv[i]);
        switch (flg){

            case str2int("-a"):
                cout << "parse : line " << __LINE__ << endl;
                ++i;
                if(i < argc) {
                    string tmp = argv[i];
                    cout << tmp << endl;
                    int sep = tmp.find(':');
                    if(sep < 0){
                        // no name is provided
                        a_usecols = explode_to_int(tmp, ',');
                        for(size_t k{}; k < a_usecols.size(); ++k){
                            a_names.push_back("<>");
                        }
                    }else {
                        string first = tmp.substr(0, sep);

                        string second = tmp.substr(sep + 1);
//                        cout << "first " << first << endl;
//                        cout << "second " << second << endl;

                        a_usecols = explode_to_int(first, ',');
                        a_names = explode_to_string(second, ',');
                        if(a_names.size() != a_usecols.size()){
                            cout << "not enough column index or name" ;
                            cout << ": line " << __LINE__ ; // comment this on deployment
                            cout << endl;
                        }
                    }
                }
                ++i;
                break;
            case str2int("-b"):
                cout << "parse : line " << __LINE__ << endl;
                ++i;
                if(i < argc) {
                    string tmp = argv[i];
                    cout << tmp << endl;
                    int sep = tmp.find(':');
                    if(sep < 0){
                        // no name is provided
                        b_usecols = explode_to_int(tmp, ',');
                        for(size_t k{}; k < b_usecols.size(); ++k){
                            b_names.push_back("<>");
                        }
                    }else {
                        string first = tmp.substr(0, sep);

                        string second = tmp.substr(sep + 1);
//                        cout << "first " << first << endl;
//                        cout << "second " << second << endl;

                        b_usecols = explode_to_int(first, ',');
                        b_names = explode_to_string(second, ',');
                        if(b_names.size() != b_usecols.size()){
                            cout << "not enough column index or name" ;
                            cout << ": line " << __LINE__ ; // comment this on deployment
                            cout << endl;
                        }
                    }
                }
                ++i;
                break;
            case str2int("-c"):
                write_header_and_comment = true;
                ++i;
                break;
            case str2int("-d"):
                ++i;
                if(i < argc) {
                    delimeter = argv[i][0];
                    cout << "delimeter provided " << int(delimeter) << endl;
//                    cout << "valid delimiter list " << endl;
//                    cout << "delimiter space = int() " << int(' ') << endl;
//                    cout << "delimiter horizontal tab = int(\\t) " << int('\t') << endl;
//                    cout << "delimiter vertical tab = int(\\v) " << int('\v') << endl;
//                    cout << "delimiter comma= int(,) " << int(',') << endl;
                }
                ++i;
                break;
            case str2int("-f"):
                cout << "parse : line " << __LINE__ << endl;
                ++i;
                if(i < argc) {
                    in_filename = argv[i];
                    cout << in_filename << endl;
                }
                ++i;
                break;
            case str2int("-i"):
                ++i;
                if(i < argc) {
                    info = argv[i];
                }
                ++i;
                break;

            case str2int("-o"):
                ++i;
                if(i < argc) {
                    out_filename = argv[i];
                }
                ++i;
                break;
            case str2int("-p"):
            case str2int("--precision"):
                ++i;
                if(i < argc) {
                    f_precision = stoi(argv[i]);
                }
                ++i;
                break;
            case str2int("-s"):
                ++i;
                if(i < argc) {
                    skiprows = stoi(argv[i]);
                }
                ++i;
                break;
            case str2int("-t"):
                ++i;
                if(i < argc) {
                    test_size = (size_t)stoi(argv[i]);
                }
                ++i;
                break;
            case str2int("--threads"):
                ++i;
                if(i < argc) {
                    n_threads = stoi(argv[i]);
                }
                ++i;
                break;
            case str2int("-h"):
            case str2int("--help"):
                help();
                exit(0);

            case str2int("-v"):
            case str2int("--version"):
                version();
                exit(0);

            case str2int("-w"):
                write_input_data = true;
                break;
            default:
                help();
                exit(0);
        }

    }

//    print_vector(a_usecols);
//    print_vector(a_names);
//    print_vector(b_usecols);
//    print_vector(b_names);

    if(out_filename.empty()){
        out_filename = in_filename + out_file_flag;
    }

    // reading input file
    if(b_usecols.empty()){
        cerr << "no specified columns" << endl;
        exit(1);
    }
//    view(b_usecols);

    delimeter = analyze_delimeter(in_filename, skiprows, delimeter);
    vector<vector<double>> b_data_in = loadtxt_v2(in_filename, b_usecols, skiprows, delimeter);
    vector<vector<double>> a_data;
    if(a_usecols.empty()){
        unsigned long N = b_data_in.size();
        unsigned long m = b_data_in[0].size();
        cout << "initializing independent data with following shape ("
             << N << "," << m << ")" << endl;
        a_data.resize(N);
        for(size_t i{}; i < N; ++i) {
            a_data[i].resize(m);
            for (size_t j{}; j < m; ++j) {
                a_data[i][j] = double (i) / N;
            }
        }
    }else {
        a_data = loadtxt_v2(in_filename, a_usecols, skiprows, delimeter);
    }
//    view_matrix(b_data_in);
    // performing convolution
    Convolution conv(n_threads);
    vector<vector<double>> b_data_out = conv.run_multi_omp(b_data_in);
//    vector<vector<double>> b_data_out = conv.run_multi_omp_v2(b_data_in);
    conv.timeElapsed();

    // writing output to file
    savetxt_multi(in_filename,
                  out_filename,
                  info,
                  write_header_and_comment,
                  delimeter,
                  write_input_data,
                  a_data,
                  b_data_in,
                  b_data_out,
                  f_precision);

}


int cmd_args_v2(int argc, char** argv){
    string in_filename;
    string out_filename;
    string out_file_flag = "_convoluted.txt";
    int header_line{0};
    size_t test_size{0};
    vector<int> a_usecols, b_usecols; // b_usecols will be convolved and a_usecols will remain unchanged
    vector<string> a_names, b_names;
    string info;
    bool write_header_and_comment{false};
    int skiprows{0};
    char delimeter{' '};
    bool write_input_data {false};
    int flg;
    int f_precision{10};
    int n_threads{1};
    double threshold{};
    int times{1};

    try
    {
        /** Define and parse the program options
         */
        namespace po = boost::program_options;
        po::options_description desc("Options");
        desc.add_options()
                ("help,h", "Print help messages")
                ("without,a", po::value<vector<int>>(&a_usecols)->multitoken()->composing(), "columns that we want in the output file without performing convolution.\nNo default value.")
                ("with,b",po::value<vector<int>>(&b_usecols)->required()->multitoken()->composing(),"columns that we want in the output file with performing convolution.")
                ("copy,c",po::value<bool>(&write_header_and_comment)->default_value(true),
                        "If provided the header and comment from the input file will be written "
                        "without modification to the output file. Header is the first line of the input file.")
                ("delimiter,d",po::value<char>()->default_value(' '), "Delimiter to use. Default value is ' '.")
                ("in",po::value<string>(&in_filename)->required(), "name of the input file that we want to convolve")
                ("out",po::value<string>(&out_filename), "name of the output file. "
                        "If not provided the string \'_convoluted.txt\' "
                        "will be appended to the input file.")
                ("info,i",po::value<string>(&info),"Info to write as comment in the output file")
                ("precision,p",po::value<int>(&f_precision)->default_value(10), "Info to write as comment in the output file")
                ("threads,t",po::value<int>(&n_threads)->default_value(1), "Info to write as comment in the output file")
                ("version,v",po::value<int>()->notifier(&version_notified), "Info to write as comment in the output file")
                ("skip",po::value<int>(&skiprows)->default_value(0), "Number of rows to skip from the input file. Default value is 0.")
                ("write,w",po::value<bool>(&write_input_data)->default_value(false),"If provided input b data will be written to the output file.")
                ("threshold",po::value<double>(&threshold)->default_value(1e-15),"If weight factor that multiplies input data at each iteration is less than\n"
                        " `threshold` then break that loop. Program performs way faster in this way.")
                ("times",po::value<int>(&times)->default_value(1),"Number of times to perform convolution.");

//        cout << __LINE__ << endl;
        po::variables_map vm;
        try   {
            po::store(po::parse_command_line(argc, argv, desc), vm); // can throw
            // there are any problems
            po::notify(vm); // throws on error, so do after help in case
            /** --help option
             */
            if (vm.count("help")) {
                std::cout << "Convolution app" << std::endl
                          << desc << std::endl;
                cout << "A line that begins with '#' is considered a commented line." << endl;
//                return SUCCESS;
            }

            cout << __LINE__ << endl;

        }
        catch(po::error& e)
        {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
            std::cerr << desc << std::endl;
            cout << __LINE__ << endl;
            return ERROR_IN_COMMAND_LINE;
        }
        catch (...) { cout << "default exception"; }
        // application code here //

    }
    catch(std::exception& e)
    {
        std::cerr << "Unhandled Exception reached the top of main: "
                  << e.what() << ", application will now exit" << std::endl;
        return ERROR_UNHANDLED_EXCEPTION;

    }
    catch (...) { cout << "default exception"; }

    if(out_filename.empty()){
        out_filename = in_filename + out_file_flag;
    }
    /*******
     * checking provided arguments
     * *****/

    cout << __LINE__ << endl;
    cout << "in_filename " << in_filename << endl;
    cout << "out_filename " << out_filename << endl;
    cout << "out_file_flag " << out_file_flag << endl;
    cout << "header_line " << header_line << endl;
    cout << "test_size " << test_size << endl;
    print_vector(a_usecols);
    print_vector(a_names);
    print_vector(b_usecols);
    print_vector(b_names);
    cout << "info " << info << endl;
    cout << "write_header_and_comment " << write_header_and_comment << endl;
    cout << "skiprows " << skiprows << endl;
    cout << "delimiter " << delimeter << endl;
    cout << "write_input_data " << write_input_data << endl;
    cout << "flg " <<  flg << endl;
    cout << "f_precision " << f_precision << endl;
    cout << "n_threads " << n_threads << endl;
    cout << "threshold " << threshold << endl;
    cout << "times " << times << endl;

    cout << __LINE__ << endl;
    delimeter = analyze_delimeter(in_filename, skiprows, delimeter);
    vector<vector<double>> b_data_in = loadtxt_v2(in_filename, b_usecols, skiprows, delimeter);
    vector<vector<double>> a_data;
    if(a_usecols.empty()){
        unsigned long N = b_data_in.size();
        unsigned long m = b_data_in[0].size();
        cout << "initializing independent data with following shape ("
             << N << "," << m << ")" << endl;
        a_data.resize(N);
        for(size_t i{}; i < N; ++i) {
            a_data[i].resize(m);
            for (size_t j{}; j < m; ++j) {
                a_data[i][j] = double (i) / N;
            }
        }
    }else {
        a_data = loadtxt_v2(in_filename, a_usecols, skiprows, delimeter);
    }
//    view_matrix(b_data_in);
    // performing convolution once
//    vector<vector<double>> b_data_out = convolve_2d(b_data_in, n_threads);
    vector<vector<double>> b_data_out = convolve_2d_fast(b_data_in, n_threads, threshold);


    // for multiple convolution
//    auto tmp = b_data_in;
//    vector<vector<double>> b_data_out;
//    for(int i{}; i < times; ++i){
//        b_data_out = convolve_2d_fast(tmp, n_threads, threshold);
//        tmp = b_data_out;
//    }

    // writing output to file
    savetxt_multi(in_filename,
                  out_filename,
                  info,
                  write_header_and_comment,
                  delimeter,
                  write_input_data,
                  a_data,
                  b_data_in,
                  b_data_out,
                  f_precision);
    return 0;
}

void help_v3(){
    string hlp = R"***(Usage:
convolution [--in <STRING>] [-a <INT>,<INT>,...] [-b <INT>,<INT>,...] [-h] [-t <INT>] [-i <STRING>]

perform convolution based on provided options.

Options                      Description
  -a, --without              columns that we want in the output file without performing convolution.
                             No default value.
  -b, --with                 columns that we want in the output file with performing convolution.
                             No default value.
  -c, --copy                 If provided the header and comment from the input file will be written
                             without modification to the output file. Header is the first line of the
                             input file.
  -d, --delimiter            Delimiter to use. Default value is ' '.
      --in                   name of the input file that we want to convolute. No default value.
  -i  --info                 Info to write as comment in the output file
      --out                  name of the output file. If not provided the string '_convoluted.txt' will be
                             appended to the input file.
  -p, --precision            Floating point precision when writing in the data file. Default value is 10
  -s, --skip                 Number of rows to skip from the input file. Default value is 0.
  -t, --threads              Explicitly specify number of thread to use. Default is the max number of thread
                             allowed by the system.
      --threshold            If weight factor that multiplies input data at each iteration is less than
                            `threshold` then break that loop. Program performs way faster in this way.
      --times                Number of times to perform convolution.
  -h, --help                 display this help and exit
  -v, --version              output version information and exit
  -w, --write                If provided input b data will be written to the output file.


The INT argument is an integer.
The STRING argument is a string of characters.

A line that begins with '#' is considered a commented line.

Exit status:
 0  if OK,
 1  if minor problems (e.g., cannot access subdirectory),
 2  if serious trouble (e.g., cannot access command-line argument).
)***";
    cout << hlp << endl;
}

int cmd_args_v3(int argc, char** argv){
    string in_filename;
    string out_filename;
    string out_file_flag = "_convoluted";
    int header_line{0};
    size_t test_size{0};
    vector<int> a_usecols, b_usecols; // b_usecols will be convolved and a_usecols will remain unchanged
    vector<string> a_names, b_names;
    string info;
    bool write_header_and_comment{false};
    int skiprows{0};
    char delimiter=' ';
    bool write_input_data {false};

    int f_precision{10};
    int n_threads{1};
    double threshold{};
    int times{1};

    parse_cmd_arg(argc, argv, in_filename, out_filename, a_usecols, b_usecols, info,
                        write_header_and_comment, skiprows, write_input_data, f_precision, n_threads,
                        threshold, times, delimiter);

    parse_cmd_arg_boost(argc, argv, in_filename, out_filename, a_usecols, b_usecols, info,
//                        write_header_and_comment, skiprows, write_input_data, f_precision, n_threads,
//                        threshold, times, delimiter);

    if(out_filename.empty()){
//        out_filename = in_filename + out_file_flag + ".txt";
        out_filename = in_filename + out_file_flag + "_" + to_string(times) + "times.txt";
    }
    /*******
     * checking provided arguments
     * *****/

    cout << __LINE__ << endl;
    cout << "in_filename " << in_filename << endl;
    cout << "out_filename " << out_filename << endl;
    cout << "out_file_flag " << out_file_flag << endl;
    cout << "header_line " << header_line << endl;
    cout << "test_size " << test_size << endl;
    print_vector(a_usecols);
    print_vector(a_names);
    print_vector(b_usecols);
    print_vector(b_names);
    cout << "info " << info << endl;
    cout << "write_header_and_comment " << write_header_and_comment << endl;
    cout << "skiprows " << skiprows << endl;
    cout << "delimiter " << delimiter << endl;
    cout << "write_input_data " << write_input_data << endl;
    cout << "f_precision " << f_precision << endl;
    cout << "n_threads " << n_threads << endl;
    cout << "threshold " << threshold << endl;
    cout << "times " << times << endl;

    cout << __LINE__ << endl;
    delimiter = analyze_delimeter(in_filename, skiprows, delimiter);
    vector<vector<double>> b_data_in = loadtxt_v2(in_filename, b_usecols, skiprows, delimiter);
    vector<vector<double>> a_data;
    if(a_usecols.empty()){
        unsigned long N = b_data_in.size();
        unsigned long m = b_data_in[0].size();
        cout << "initializing independent data with following shape ("
             << N << "," << m << ")" << endl;
        a_data.resize(N);
        for(size_t i{}; i < N; ++i) {
            a_data[i].resize(m);
            for (size_t j{}; j < m; ++j) {
                a_data[i][j] = double (i) / N;
            }
        }
    }else {
        a_data = loadtxt_v2(in_filename, a_usecols, skiprows, delimiter);
    }
//    view_matrix(b_data_in);

    // performing convolution once
//    vector<vector<double>> b_data_out = convolve_2d(b_data_in, n_threads);
//    vector<vector<double>> b_data_out = convolve_2d_fast(b_data_in, n_threads, threshold);


    // for multiple convolution
    auto tmp = b_data_in;
    vector<vector<double>> b_data_out;
    for(int i{}; i < times; ++i){
        if(threshold == -1) {
            b_data_out = convolve_2d(b_data_in, n_threads);
        }else {
            b_data_out = convolve_2d_fast(tmp, n_threads, threshold);
        }
        tmp = b_data_out;
    }

    // writing output to file
    savetxt_multi(in_filename,
                  out_filename,
                  info,
                  write_header_and_comment,
                  delimiter,
                  write_input_data,
                  a_data,
                  b_data_in,
                  b_data_out,
                  f_precision);
    return 0;
}

int parse_cmd_arg_boost(int argc, char *const *argv, string &in_filename, string &out_filename, vector<int> &a_usecols,
                         vector<int> &b_usecols, string &info, bool &write_header_and_comment, int &skiprows,
                         bool &write_input_data, int &f_precision, int &n_threads, double &threshold, int &times, char& delimiter) {

    write_input_data=false;
    write_header_and_comment = true;

    try
    {
        /** Define and parse the program options
         */
        namespace po = boost::program_options;
        boost::program_options::options_description desc("Options");
        desc.add_options()
                ("help,h", "Print help messages")
                ("without,a", boost::program_options::value<vector<int>>(&a_usecols)->multitoken()->composing(), "columns that we want in the output file without performing convolution.\nNo default value.")
                ("with,b", boost::program_options::value<vector<int>>(&b_usecols)->required()->multitoken()->composing(), "columns that we want in the output file with performing convolution.")
                ("copy,c", "If provided the header and comment from the input file will be written "
                         "without modification to the output file. Header is the first line of the input file.")
                ("delimiter,d", boost::program_options::value<char>(&delimiter)->default_value(' '), "Delimiter to use. Default value is ' '.")
                ("in", boost::program_options::value<string>(&in_filename)->required(), "name of the input file that we want to convolve")
                ("out", boost::program_options::value<string>(&out_filename), "name of the output file. "
                        "If not provided the string \'_convoluted.txt\' "
                        "will be appended to the input file.")
                ("info,i", boost::program_options::value<string>(&info), "Info to write as comment in the output file")
                ("precision,p", boost::program_options::value<int>(&f_precision)->default_value(10), "Info to write as comment in the output file")
                ("threads,t", boost::program_options::value<int>(&n_threads)->default_value(1), "Info to write as comment in the output file")
                ("version,v", boost::program_options::value<int>()->notifier(&version_notified), "Info to write as comment in the output file")
                ("skip", boost::program_options::value<int>(&skiprows)->default_value(0), "Number of rows to skip from the input file. Default value is 0.")
                ("write,w", "If provided input b data will be written to the output file.")
                ("threshold", boost::program_options::value<double>(&threshold)->default_value(1e-15), "If weight factor that multiplies input data at each iteration is less than\n"
                        " `threshold` then break that loop. Program performs way faster in this way.")
                ("times", boost::program_options::value<int>(&times)->default_value(1), "Number of times to perform convolution.");

//        cout << __LINE__ << endl;
        boost::program_options::variables_map vm;
        try   {
            store(parse_command_line(argc, argv, desc), vm); // can throw
            // there are any problems
            notify(vm); // throws on error, so do after help in case
            /** --help option
             */
            if (vm.count("help")) {
                cout << "Convolution app" << endl
                          << desc << endl;
                cout << "A line that begins with '#' is considered a commented line." << endl;
                return SUCCESS;
            }
            if (vm.count("write") || vm.count("w")) {
                write_input_data = true;
                cout << "input data will be written" << endl;
//                return SUCCESS;
            }
            if (vm.count("copy")|| vm.count("c")) {
                write_header_and_comment = false;
                cout << "header information will not be written" << endl;
//                return SUCCESS;
            }

            cout << __LINE__ << endl;

        }
        catch(boost::program_options::error& e)
        {
            cerr << "ERROR: " << e.what() << endl << endl;
            cerr << desc << endl;
            cout << __LINE__ << endl;
            return ERROR_IN_COMMAND_LINE;
        }
        catch (...) { cout << "default exception"; }
        // application code here //

    }
    catch(exception& e)
    {
        cerr << "Unhandled Exception reached the top of main: "
                  << e.what() << ", application will now exit" << endl;
        return ERROR_UNHANDLED_EXCEPTION;

    }
    catch (...) { cout << "default exception"; }
}

void parse_cmd_arg(int argc, char *const *argv, string &in_filename, string &out_filename, vector<int> &a_usecols,
                         vector<int> &b_usecols, string &info, bool &write_header_and_comment, int &skiprows,
                         bool &write_input_data, int &f_precision, int &n_threads, double &threshold, int &times, char& delimiter) {
    if(argc == 1){
        help_v3();
        exit(0);
    }
    int flg;
    for(int i{1}; i < argc;){
        cout << argv[i] << endl;

        flg = str2int(argv[i]);
        switch (flg){

            case str2int("-a"):
            case str2int("--without"):
                cout << "parse : line " << __LINE__ << endl;
                ++i;
                if(i < argc) {
                    string tmp = argv[i];
                    cout << tmp << endl;
                    int sep = tmp.find(':');
                    if(sep < 0){
                        // no name is provided
                        a_usecols = explode_to_int(tmp, ',');
                    }else {
                        string first = tmp.substr(0, sep);

                        string second = tmp.substr(sep + 1);
//                        cout << "first " << first << endl;
//                        cout << "second " << second << endl;

                        a_usecols = explode_to_int(first, ',');
                    }
                }
                ++i;
                break;
            case str2int("-b"):
            case str2int("--with"):
                cout << "parse : line " << __LINE__ << endl;
                ++i;
                if(i < argc) {
                    string tmp = argv[i];
                    cout << tmp << endl;
                    int sep = tmp.find(':');
                    if(sep < 0){
                        // no name is provided
                        b_usecols = explode_to_int(tmp, ',');

                    }else {
                        string first = tmp.substr(0, sep);

                        string second = tmp.substr(sep + 1);
//                        cout << "first " << first << endl;
//                        cout << "second " << second << endl;

                        b_usecols = explode_to_int(first, ',');

                    }
                }
                ++i;
                break;
            case str2int("-c"):
            case str2int("--copy"):
                write_header_and_comment = true;
                ++i;
                break;
            case str2int("-d"):
            case str2int("--delimiter"):
                ++i;
                if(i < argc) {
                    delimiter = argv[i][0];
                    cout << "delimeter provided " << int(delimiter) << endl;
//                    cout << "valid delimiter list " << endl;
//                    cout << "delimiter space = int() " << int(' ') << endl;
//                    cout << "delimiter horizontal tab = int(\\t) " << int('\t') << endl;
//                    cout << "delimiter vertical tab = int(\\v) " << int('\v') << endl;
//                    cout << "delimiter comma= int(,) " << int(',') << endl;
                }
                ++i;
                break;
            case str2int("--in"):
                cout << "parse : line " << __LINE__ << endl;
                ++i;
                if(i < argc) {
                    in_filename = argv[i];
                    cout << in_filename << endl;
                }
                ++i;
                break;
            case str2int("-i"):
            case str2int("--info"):
                ++i;
                if(i < argc) {
                    info = argv[i];
                }
                ++i;
                break;

            case str2int("-o"):
            case str2int("--out"):
                ++i;
                if(i < argc) {
                    out_filename = argv[i];
                }
                ++i;
                break;
            case str2int("-p"):
            case str2int("--precision"):
                ++i;
                if(i < argc) {
                    f_precision = stoi(argv[i]);
                }
                ++i;
                break;
            case str2int("--skip"):
                ++i;
                if(i < argc) {
                    skiprows = stoi(argv[i]);
                }
                ++i;
                break;
            case str2int("-t"):
            case str2int("--threads"):
                ++i;
                if(i < argc) {
                    n_threads = stoi(argv[i]);
                }
                ++i;
                break;
            case str2int("-h"):
            case str2int("--help"):
                help_v3();
                exit(0);

            case str2int("-v"):
            case str2int("--version"):
                version();
                exit(0);

            case str2int("-w"):
            case str2int("--write"):
                write_input_data = true;
                break;
            case str2int("--times"):
                ++i;
                times = stoi(argv[i]);
                ++i;
                break;
            default:
                help_v3();
                exit(0);
        }

    }
}