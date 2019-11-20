//
// Created by shahnoor on 11/20/19.
//

#include "process.h"
#include "../include/string_methods.h"
#include <iostream>

using namespace std;


void test_process(int argc, char **argv){
    test_process_v2(argc, argv);
}

void test_process_v1(int argc, char **argv){
    for(int i{}; i < argc; ++i){
        cout << argv[i] << endl;
        if(int(argv[i][0]) == int('-')){ // options
            ++i;
            cout << "while" << endl;
            while(i < argc && int(argv[i][0]) != int('-')){

                cout << argv[i] << endl;

                ++i;
            }
            --i;
        }
    }
}

void test_process_v2(int argc, char **argv){
    for(int i{}; i < argc; ++i){
        cout << argv[i] << endl;
        if(int(argv[i][0]) == int('-')){ // options
            ++i;
            auto arr = get_int_array(argc, argv, i);
            cout << "{";
            for(auto a: arr){
                cout << a << ",";
            }
            cout << "}" << endl;
        }
    }
}


bool is_number(string str){
//    cout << "checking " << str << " size " << str.size() << endl;
    for(int i = 0;i < str.size();i++) {
//        cout  << str[i] << endl;
        if(!isdigit(str[i])) {
            return false;
        }
    }
    return true;
}

//bool is_number(char* str){
//    cout << "checking [" << str << "] size " << strlen(str) << endl;
//    for(int i = 0;i < strlen(str);i++) {
//        cout  << str[i] << endl;
//        if(!isdigit(str[i])) {
//            return false;
//        }
//    }
//    return true;
//}



vector<int> get_int_array(int argc, char **argv, int &index){
    vector<int> arr;
    cout << "while" << endl;
    while(index < argc && int(argv[index][0]) != int('-')){

        cout << argv[index] << endl;
        if(!is_number(argv[index])){
            cout << "NAN error" << endl;
            break;
        }else {
            auto num = stoi(argv[index]);
            arr.emplace_back(num);
        }
        ++index;
    }
    --index;

    return arr;
}