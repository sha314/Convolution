//
// Created by shahnoor on 9/14/18.
//

#include <algorithm>
#include "include/string_methods.h"

std::string trim(std::string str, char ch) {
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), [&](int c) {
        return c != ch;
    }));
    return str;
}


/**
 * string to int conversion. using loop.
 * @param str
 * @return
 */
unsigned int str_to_int(const char* str)
{
    unsigned x = 1;
    int len = strlen(str);
    // cout << "len, str " << len << ", " << str << endl;
    for(int i=len-1; i>=0 ; --i){
        x *= 2;
        x = x ^ str[i];
    }
    // cout << x << endl;
    return x;
}