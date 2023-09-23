#include "py_bri_wrapper.h"
#include "bri_get.h"

static PyObject* py_initialize_bri(PyObject *self, PyObject *args)
{
    const char *bam_file;
    const char *index_file;
    
    if (!PyArg_ParseTuple(args, "ss", &bam_file, &index_file)) {
        return NULL;
    }

    bam_read_idx *bri = initialize_bri(bam_file, index_file);
    if (!bri) {
        PyErr_SetString(PyExc_RuntimeError, "Could not initialize bri.");
        return NULL;
    }

    // Assuming bam_read_idx is fully encapsulated in C and not exposed to Python
    return Py_BuildValue("K", (unsigned long long) bri);
}

static PyObject* py_query_by_readname(PyObject *self, PyObject *args)
{
    unsigned long long bri_ptr;
    PyObject *read_names_list;
    const char *input_bam;
    // const char *output_sam;
    if (!PyArg_ParseTuple(args, "ssKO", &input_bam, &bri_ptr, &read_names_list)) {
        return NULL;
    }
    
    if (!PyList_Check(read_names_list)) {
        PyErr_SetString(PyExc_TypeError, "read_names must be a list.");
        return NULL;
    }

    bam_read_idx *bri = (bam_read_idx *) bri_ptr;
    Py_ssize_t size = PyList_Size(read_names_list);
    char **read_names = malloc(size * sizeof(char *));
    
    for (Py_ssize_t i = 0; i < size; i++) {
        PyObject *item = PyList_GetItem(read_names_list, i);
        if (!PyUnicode_Check(item)) {
            free(read_names);
            PyErr_SetString(PyExc_TypeError, "All read names must be strings.");
            return NULL;
        }
        
        read_names[i] = PyUnicode_AsUTF8(item);
    }

    // query_by_readname(input_bam, bri, read_names, (uint64_t)size, output_sam);
    kstring_t ks = {0, 0, NULL};  // 初始化 kstring_t

    query_by_readname(input_bam, bri, read_names, (uint64_t)size, &ks);

    PyObject *result = Py_BuildValue("s", ks.s);  // 将 kstring_t 转换为 Python 字符串

    free(ks.s);  // 释放 kstring_t
    free(read_names);  // 释放 read_names

    return result;  // 返回 Python 字符串
}

static PyMethodDef BriMethods[] = {
    {"initialize_bri",  py_initialize_bri, METH_VARARGS, "Initialize bam_read_idx."},
    {"query_by_readname",  py_query_by_readname, METH_VARARGS, "Query by read name."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef brimodule = {
    PyModuleDef_HEAD_INIT,
    "bri",
    NULL,
    -1,
    BriMethods
};

PyMODINIT_FUNC PyInit_bri(void)
{
    return PyModule_Create(&brimodule);
}
