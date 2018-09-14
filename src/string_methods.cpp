//
// Created by shahnoor on 9/14/18.
//

#include <algorithm>
#include "string_methods.h"

std::string trim(std::string str, char ch) {
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), [&](int c) {
        return c != ch;
    }));
    return str;
}
