//
// Created by shahnoor on 11/19/19.
//

#ifndef CONVOLUTION_LOGGER_H
#define CONVOLUTION_LOGGER_H

/**
 * Writes log to a file
 */
#include <fstream>
#include <utility>
#include <vector>

/**
 Only one instance of Logging class is allowed. No new instance can be created.
 Whenever it is called all the previous changes will be present.
 Use :
 	Database integration.
 	Saving configuration or setting.
 	Logging data.
*/
class Logging{
    static Logging* _instance;
    std::string _filename;
    std::ofstream _fout;

    Logging(){
        _filename="a.txt";
        _fout.open(_filename);
    }

public:
    static Logging* getInstance(){
        if (_instance == nullptr){
            _instance = new Logging();
        }
        return _instance;
    }

    void setFilename(std::string filename){
        _filename= std::move(filename);
        _fout.close();
        _fout.open(_filename);
    }
    void addText(const std::string &text){
        _fout << text << std::endl;
    }

    void addText(double text){
        _fout << text << std::endl;
    }

    void addText(std::vector<double>& text){
        for(auto a: text) {
            _fout << a << "\t";
        }
        _fout << std::endl;
    }

};

#endif //CONVOLUTION_LOGGER_H
