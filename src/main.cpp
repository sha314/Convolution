#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <map>
#include <chrono>
#include <iomanip>


#include "convolution/binomial.h"
#include "convolution/convolution.h"
#include "io/data_reader.h"
#include "include/printer.h"
#include "include/string_methods.h"
#include "io/data_writer.h"
#include "tests/test1.h"
#include "tests/test2.h"
#include "cmd_args.h"
#include "io/logger.h"
#include "args/process.h"


using namespace std;



/**
 * Can store any type of value and show it later.
 * But how to return it
 */
class MultiType{
    int val_int;
    double val_double;
    string val_string;
    vector<int> arr_int;
    vector<double> arr_double;
    vector<string> arr_string;

    vector<int> value_provided={-1,-2,-3,-4,-5,-6};
public:
    ~MultiType() = default;
    MultiType(int a){
        val_int = a;
        value_provided[0] = 1;
    }
    MultiType(double a){
        val_double = a;
        value_provided[1] = 1;
    }
    MultiType(string a){
        val_string = a;
        value_provided[2] = 1;
    }
    MultiType(vector<int>  a){
        arr_int = std::move(a);
        value_provided[3] = 1;
    }
    MultiType(vector<double>  a){
        arr_double = std::move(a);
        value_provided[4] = 1;
    }
    MultiType(vector<string> a){
        arr_string = std::move(a);
        value_provided[5] = 1;
    }

    void show(){
        if(value_provided[0] == 1) cout << val_int;
        if(value_provided[1] == 1) cout << val_double;
        if(value_provided[2] == 1) cout << val_string;
        if(value_provided[3] == 1) {
            cout << "{";
            for(size_t i{}; i < arr_int.size(); ++i){
                cout << arr_int[i] <<",";
            }
            cout << "}";
        }
        if(value_provided[4] == 1) {
            cout << "{";
            for(size_t i{}; i < arr_double.size(); ++i){
                cout << arr_double[i] <<",";
            }
            cout << "}";
        }
        if(value_provided[5] == 1) {
            cout << "{";
            for(size_t i{}; i < arr_string.size(); ++i){
                cout << arr_string[i] <<",";
            }
            cout << "}";
        }
    }



};

/**
 * Takes command line arguments as constructing
 * Processes
 */
class ArgOperator{

public:
    ~ArgOperator() = default;
    ArgOperator(int argc, char* argv[]);
    void parse_arguments();

};


void on_age(int age)
{
    std::cout << "On age: " << age << '\n';
}


Logging* Logging::_instance=nullptr;

/***
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
    cout << "Convolution of big data" << endl;

    auto t0 = std::chrono::system_clock::now();

//    cmd_args(argc, argv);
//    cmd_args_v2(argc, argv);
//    cmd_args_v3(argc, argv);
//    test1_convolution();
//    test2_convolution();
//      test3_convolution();
//    test4_convolution();
    test_process(argc, argv);

    auto t1 = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = t1-t0;
    std::time_t end_time = std::chrono::system_clock::to_time_t(t1);
#ifdef DEBUG_FLAG
    cout << "flag DEBUG_FLAG is defined" << endl;
#endif
#ifdef UNIT_TEST
    cout << "flag UNIT_TEST is defined" << endl;
#endif
#ifdef USE_BOOST
    // if provided during compile time boost part will be selected for compilation
    cout << "flag USE_BOOST is defined" << endl;
#endif
    cout << "Program finished at " << std::ctime(&end_time) << endl;
    cout << "Total Time elapsed " << elapsed_seconds.count()/60 << " minutes" << endl;
    return 0;
}

