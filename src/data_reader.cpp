//
// Created by shahnoor on 2/2/2018.
//

#include "data_reader.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>


using namespace std;



/**
 * Requirement:
 *  Only First line of the file should contain the header information
 *  // todo start scanning from '{' and end at '}'
 * @param filename
 * @param delimiter
 * @param comment
 * @return
 */
std::map<std::string, unsigned>
read_header_json(std::string filename, char comment)
{

    ifstream fin(filename);

    map<string, unsigned> header_info;
    string line;
    string key;
    unsigned value;
    bool flag{false};
    while (getline(fin, line)){
        if(line[0] == comment) {
            continue;
        }else{
            cout << "Header -> ";
            cout << line << endl;
            break;
        }
    }

    // removing all spaces
    line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

    size_t a = line.find('{')+1;
    size_t b = line.find('}')-1;
//    cout << a << " to " << b << endl;
    string contents = line.substr(a, b);
//    cout << contents << endl;
    vector<string> pairs = explode_to_string(contents, ',');
    for(string &s : pairs){
        vector<string> part = explode_to_string(s, ':');

        string key_tmp = explode_to_string(part[0], '\"')[0];
        unsigned value_tmp = atoi(part[1].c_str());

        header_info[key_tmp] = value_tmp;

    }


    return header_info;
}

/**
 * Divide a string to array of string with respect to delimeter
 * @param str
 * @param ch
 * @return
 */
std::vector<std::string> explode_to_string(const std::string &str, const char &ch) {
    std::string next;
    vector<std::string> result;

    // For each character in the string
    for (std::string::const_iterator it = str.begin(); it != str.end(); it++) {
        // If we've hit the terminal character
        if (*it == ch) {
            // If we have some characters accumulated
            if (!next.empty()) {
                // Add them to the result vector
                result.push_back(next);
                next.clear();
            }
        } else {
            // Accumulate the next character into the sequence
            next += *it;
        }
    }
    if (!next.empty())
        result.push_back(next);
    return result;
}



/**
 * explode a string to vector of string
 * @param s
 * @param c
 * @return
 */
vector<int> explode_to_int(const string &s, const char &c)
{
    string buff{""};
    vector<int> v;

    for(auto n:s)
    {
        if(n != c) buff+=n; else
        if(n == c && buff != "") { v.push_back(atoi(buff.c_str())); buff = ""; }
    }
    if(buff != "") v.push_back(atoi(buff.c_str()));

    return v;
}

/**
 * Reads data from files
 * @param filename  : name of the file
 * @param usecols   : column to read
 * @param skiprows  : number of rows to be skipped (commented or uncommented)
 * @param delimiter : character used as delemeter in the file
 * @param comment   : character used as comment in the file
 * @return : data of the column
 */
vector<double> loadtxt(string filename, int usecol,
                       int skiprows, char delimiter, char comment){
    vector<double> data;
    ifstream fin(filename);

    string line;
    double value;
    unsigned r{}, c{};
    while (getline(fin, line)){
        if(r < skiprows){
            ++r;
            continue;
        }
        if(line[0] == comment) {
            continue;
        }
        istringstream iss(line);
        c = 0;
        while(iss >> value) {
            if(c == usecol){
                break;
            }
        }

        data.push_back(value);
//        cout << line << endl;
    }

    return data;
}



/**
 * Reads columns of data from files
 * @param filename  : name of the file
 * @param usecols   : columns to read
 * @param skiprows  : number of rows to be skipped (commented or uncommented)
 * @param delimiter : character used as delemeter in the file
 * @param comment   : character used as comment in the file
 * @return : data of the columns
 */
vector<vector<double>> loadtxt(string filename, const vector<int>& usecols,
                               int skiprows, char delimiter, char comment){
    vector<vector<double>> data;
    ifstream fin(filename);

    vector<double> tmp, filtered;
    string line;
    double value;
    unsigned r{};
    while (getline(fin, line)){
        if(r < skiprows){
            ++r;
            continue;
        }
        if(line[0] == comment) {
            continue;
        }
        istringstream iss(line);
        while(iss >> value) {
            tmp.push_back(value);
//            cout << value << '\t';
        }
//        cout << endl;
        for(auto c : usecols){
            filtered.push_back(tmp[c]);
        }
        data.push_back(filtered);
        tmp.clear();
        filtered.clear();
//        cout << line << endl;
    }

    return data;
}



/**
 * Get a header for output file in raw format
 * @param icolumn_name
 * @param usecols_names
 * @param header
 * @return
 */
string
output_header_raw(
        const string& icolumn_name,
        const vector<string>& usecols_names,
        map<string, unsigned>& header
){
    ostringstream oss;

    oss << '#' << icolumn_name;
    for(size_t i{0}; i != usecols_names.size(); ++i){
        oss << "\t<" + usecols_names[i] + ">";
        oss << "\t<" + usecols_names[i] + " convolved>";
    }

    oss << "#Convolved data" << endl;
    oss << "BEGIN_HEADER" << endl;
    oss << "ensemble_size\t" << header["ensemble_size"] << endl;
    oss << "length\t" << header["length"] << endl;
    oss << "data_line\t" << 8 << endl;
    oss << "END_HEADER" << endl;
    return oss.str();
}

/**
 * Get a header for output file in JSON format
 * @param icolumn_name
 * @param usecols_names
 * @param header
 * @return
 */
string output_header_json(
        const string& icolumn_name,
        const vector<string>& usecols_names,
        const map<string, unsigned>& header
){
    ostringstream oss;
    size_t max = header.size() - 1;
    int i{};
    oss << '{';
    for(auto x: header){
        oss << "\"" << x.first << "\"" << ':' << x.second ;
        if (i < max){
            oss<< ',';
        }
        ++i;
    }
    oss << '}' << endl;

    oss << "#Convolved data" << endl;

    oss << '#' << icolumn_name;
    for(size_t i{0}; i != usecols_names.size(); ++i){
        oss << "\t<" + usecols_names[i] + ">";
        oss << "\t<" + usecols_names[i] + " convolved>";
    }

    return oss.str();
}


