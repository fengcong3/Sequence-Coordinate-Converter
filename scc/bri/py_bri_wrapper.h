#ifndef PY_BRI_WRAPPER_H
#define PY_BRI_WRAPPER_H

#include <Python.h>
#include "bri_get.h"

// 声明日志函数
void debug_log(const char* format, ...);

// 声明主要的Python接口函数
static PyObject* py_initialize_bri(PyObject *self, PyObject *args);
static PyObject* py_create_server(PyObject *self, PyObject *args);
static PyObject* py_connect_to_server(PyObject *self, PyObject *args);
static PyObject* py_disconnect_from_server(PyObject *self, PyObject *args);
static PyObject* py_stop_server(PyObject *self, PyObject *args);
static PyObject* py_query_by_readname(PyObject *self, PyObject *args);

// 声明内部辅助函数
size_t calculate_total_size(bam_read_idx *bri);
void* setup_direct_access_mode(void *shm_addr);
void validate_shm_structures(void *shm_addr);

#endif // PY_BRI_WRAPPER_H
