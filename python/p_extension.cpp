//
// Created by shahnoor on 8/6/18.
//


#include "./include/p_extension.h"
#include "./include/converter.h"
#include "../src/include/convolution.h"

// use with python3.x
using namespace std;

struct module_state {
    PyObject *error;
};


#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

void view(vector<vector<double>>& mat){
    cout << "view matrix in c++" << endl;
    cout << '{';
    for(size_t i{}; i < mat.size(); ++i){
        cout << '{';
        for(size_t j{}; j < mat[i].size(); ++j){
            cout << mat[i][j] << ',';
        }
        cout << '}' << endl;
    }
    cout << '}' << endl;
    cout << "exiting C++" << endl;
}

/**
copy imput python list and convert it to c++ vector
then square each element it and return the squared sum
*/
PyObject *
view_matrix(PyObject *self, PyObject *list)
{
    cout << "entered C++" << endl;
    vector<vector<double>> mat;
    double value;
    if(PyList_Check(list)){
        Py_ssize_t i,j, r, c;
        PyObject *item, *next_item;
        r = PyList_Size(list);
        mat.resize(r);
        // cout << "number of elements " << n << endl;
        if (r < 0)  return nullptr; /* Not a list */
        for (i = 0; i < r; i++) {
            //      cout << "line :" << __LINE__ << endl;
            item = PyList_GetItem(list, i); /* Can't fail */
            if(PyList_Check(list)){
                c = PyList_Size(item);
                if (c < 0) return nullptr; /* Not a list */
                mat[i].resize(c);
                for (j = 0; j < c; j++) {
                    next_item = PyList_GetItem(item, j); /* Can't fail */
                    value = PyFloat_AsDouble(next_item);
                    if (value == -1 && PyErr_Occurred()){
                        /* Integer too big to fit in a C double, bail out */
                        std::cout << "Integer too big to fit in a C double, bail out. line :" << __LINE__ << std::endl;
                        exit(-1);
                    }
                    mat[i][j] = value;
                }

            }
        }
    }
    view(mat);

    return Py_BuildValue("l", 0);
}


static PyObject*
convolute(PyObject *self, PyObject *args, PyObject *kwargs){

    static char* keywords[] = {"a", NULL}; // "a", "b" are the keyword of the argument
    PyObject* a;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", keywords, &a)) {
        return NULL;
    }
    cout << "c++ : parsing first list" << endl;
    auto mat1 = get_matrix(a);
    view(mat1);
    Convolution conv;
    cout << "executing conv" << endl;
    auto mat_out = conv.run_multi_omp(mat1);
    cout << "got mat out" << endl;
    view(mat_out);
    PyObject* result = get_python_matrix(mat_out);
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
        {"view_matrix", (PyCFunction)view_matrix, METH_O, "give a MxN matrix and I will display it"},
        {"convolute", (PyCFunction)convolute, METH_VARARGS | METH_KEYWORDS, "matrix multiplicatoin"},
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
PyInit_convolution(void)
{

    PyObject *module = PyModule_Create(&moduledef);


    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("convolution.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

    return module;

}
