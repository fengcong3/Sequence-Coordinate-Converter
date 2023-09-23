#ifndef PY_BRI_WRAPPER_H
#define PY_BRI_WRAPPER_H

#include <Python.h>

// Python method definitions
static PyObject* py_initialize_bri(PyObject *self, PyObject *args);
static PyObject* py_query_by_readname(PyObject *self, PyObject *args);

#endif // PY_BRI_WRAPPER_H
