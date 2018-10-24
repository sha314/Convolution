//
// Created by shahnoor on 8/6/18.
//

#ifndef CONVOLUTION_P_EXTENSION_H
#define CONVOLUTION_P_EXTENSION_H

#include <Python.h>
#include <iostream>
#include <vector>

// /****
// Numpy array as input
// *****/
// PyObject * view_numpy_array(PyObject *self, PyArrayObject *np_array);

void view(std::vector<std::vector<double>>& mat);

static PyObject * error_out(PyObject *m);

static PyObject* view_matrix(PyObject *self, PyObject *args);
static PyObject* convolute(
        PyObject *self, PyObject *args, PyObject *kwargs);


#endif //CONVOLUTION_P_EXTENSION_H
