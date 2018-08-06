//
// Created by shahnoor on 8/6/18.
//

#include "p_extension.h"

extern "C"{
#include <python3.6/Python.h>
}

#include <iostream>
#include "../convolution.h"

// use with python3.x
using namespace std;

struct module_state {
    PyObject *error;
};


#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))


static PyObject *
convolution_py(PyObject *self, PyObject *args) {
    PyObject *result;
    size_t n;
    if (!PyArg_ParseTuple(args, "L", &n))
        return NULL;
    cout << "Calculating Convolution Inside C++" << endl;
    size_t out = 0;
    result = Py_BuildValue("i", out);
    return result;
}

static PyObject *
error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

/**
List of method definitions
**/
static PyMethodDef myextension_methods[] = {
        {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
        {"convolution", (PyCFunction)convolution_py, METH_VARARGS, "calculates factorial"},
        {NULL, NULL}
};


static int myextension_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int myextension_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

/**
module definition
**/
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "convolution",
        NULL,
        sizeof(struct module_state),
        myextension_methods,
        NULL,
        myextension_traverse,
        myextension_clear,
        NULL
};

#define INITERROR return NULL

/**
 module initialization
**/
PyMODINIT_FUNC
PyInit_myextension(void)
{

    PyObject *module = PyModule_Create(&moduledef);


    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("myextension.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

    return module;

}