//
// Created by shahnoor on 10/27/18.
//

#ifndef CONVOLUTION_UTYPE_H
#define CONVOLUTION_UTYPE_H

#include <cstddef>
#include <string>

/***********************
 * Universal Type class
 * An object capable of holding different type of data but not at the same time.
 *
 */
class UType{
    int _d_int;
    long _d_long;
    size_t _d_size_t;
    double _d_double;
    float _d_float;

    std::string _type;
public:
    ~UType() = default;
    UType(int x);
    UType(long x);
    UType(float x);
    UType(double x);
    UType(size_t x);

    std::string type() const {return _type;}

    void get(int& x);
    void get(long& x);
    void get(float& x);
    void get(double& x);
    void get(size_t& x);

};


#endif //CONVOLUTION_UTYPE_H
