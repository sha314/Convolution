//
// Created by shahnoor on 11/19/19.
//

#ifndef CONVOLUTION_ARRAY_H
#define CONVOLUTION_ARRAY_H

#include <vector>
#include <cstddef>

namespace num_array {
    template<typename T>
    T max(const std::vector <T> &vec) {
        T a = vec[0];
        for (auto b: vec){
            if(b > a){
                a = b;
            }
        }
        return a;
    }

    template<typename T>
    T max(const std::vector <std::vector<T>> &vec, size_t column) {
        T a = vec[0][column];
        for(size_t row{}; row < vec.size(); ++row){
            auto b = vec[row][column];
            if(b > a){
                a = b;
            }
        }
        return a;
    }

    template<typename T>
    T min(const std::vector <T> &vec) {
        T a = vec[0];
        for (auto b: vec){
            if(b < a){
                a = b;
            }
        }
        return a;
    }

    template<typename T>
    size_t argmax(const std::vector <T> &vec) {
        T a = vec[0];
        size_t index=0;
        for(size_t i{}; i < vec.size(); ++i){
            auto b = vec[i];
            if(b > a){
                a = b;
                index = i;

            }
        }
        return index;
    }

    template<typename T>
    size_t argmin(const std::vector <T> &vec) {
        T a = vec[0];
        size_t index=0;
        for(size_t i{}; i < vec.size(); ++i){
            auto b = vec[i];
            if(b < a){
                a = b;
                index = i;

            }
        }
        return index;
    }
 /******************************
 * Perform differentiation
 *******************************/
/**
 * Differentiation of a single column
 * @tparam T
 * @param vec
 * @param column
 * @return
 */
    template<typename T>
    std::vector<T> diff(const std::vector <std::vector<T>> &vec, size_t column){
        T prev = vec[0][column];
        std::vector<T> arr(vec.size()); // number of rows must be equal
        for(size_t row{1}; row < vec.size(); ++row){
            T now = vec[row][column];
            arr[row] = now - prev;
            prev = now; // replace old value
        }
        return arr;
    }
}

#endif //CONVOLUTION_ARRAY_H
