//
// Created by shahnoor on 9/14/18.
//

#ifndef CONVOLUTION_STRING_METHODS_H
#define CONVOLUTION_STRING_METHODS_H

#include <string>
#include <cstring>

std::string trim(std::string str, char ch);
unsigned int str_to_int(const char* str);

/**
 * string to int conversion. using recursion.
 * @param str
 * @param h
 * @return
 */
constexpr unsigned int str2int(const char* str, int h=0)
{
    return !str[h] ? 1 : (str2int(str, h+1) * 2) ^ str[h];
}



#endif //CONVOLUTION_STRING_METHODS_H
