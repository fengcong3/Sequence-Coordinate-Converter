#include "py_bri_wrapper.h"
#include "bri_get.h"
#include <sys/shm.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>

#define SHM_KEY_PATH "/tmp/scc_shm_key"
#define SHM_SIZE_PATH "/tmp/scc_shm_size"
#define SHM_PROJ_ID 65
#define DEBUG_LOG 0  // 启用调试日志

// 增加日志函数
void debug_log(const char* format, ...) {
    if (!DEBUG_LOG) return;
    
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[BRI-DEBUG] ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
    fflush(stderr);
}

// 共享内存结构的头部
typedef struct {
    size_t total_size;          // 共享内存总大小
    size_t bri_offset;          // bri结构体的偏移量
    size_t readnames_offset;    // readnames数据的偏移量
    size_t records_offset;      // records数据的偏移量
    size_t name_count_bytes;    // readnames的大小
    size_t record_count;        // records的数量
    int magic_number;           // 魔数，用于验证共享内存结构的有效性
    char bam_path[1024];        // bam文件路径
    int new_approach;           // 标记使用新方法，避免与旧版本冲突
    int direct_access_mode;     // 标记使用直接访问模式
} shm_header_t;

#define SHM_MAGIC_NUMBER 0x42524953  // "BRIS" in ASCII

// 注意：移除这里对client_context_t的重复定义
// 使用bri_get.h中定义的client_context_t

static PyObject* py_initialize_bri(PyObject *self, PyObject *args)
{
    const char *bam_file;
    const char *index_file;
    
    if (!PyArg_ParseTuple(args, "ss", &bam_file, &index_file)) {
        return NULL;
    }

    debug_log("Initializing BRI for %s with index %s", bam_file, index_file);
    bam_read_idx *bri = initialize_bri(bam_file, index_file);
    if (!bri) {
        PyErr_SetString(PyExc_RuntimeError, "Could not initialize bri.");
        return NULL;
    }
    debug_log("BRI initialized successfully: %p", (void*)bri);
    debug_log("BRI stats: readnames=%p, records=%p, name_count=%zu, record_count=%zu", 
              bri->readnames, bri->records, bri->name_count_bytes, bri->record_count);

    return Py_BuildValue("K", (unsigned long long) bri);
}

// 获取共享内存总大小
size_t calculate_total_size(bam_read_idx *bri) {
    size_t header_size = sizeof(shm_header_t);
    size_t bri_size = sizeof(bam_read_idx);
    size_t readnames_size = bri->name_count_bytes;
    size_t records_size = bri->record_count * sizeof(bam_read_idx_record);
    
    size_t total = header_size + bri_size + readnames_size + records_size;
    debug_log("Calculated shared memory size: %zu bytes (header=%zu, bri=%zu, names=%zu, records=%zu)", 
              total, header_size, bri_size, readnames_size, records_size);
    
    return total;
}

// 使用extern声明这些函数，实现留在bri_get.c中
extern const char* get_record_readname(client_context_t *ctx, size_t index);
extern int direct_search_record_by_name(client_context_t *ctx, const char *readname, size_t *start_idx, size_t *end_idx);

// 修改重建指针函数，采用直接访问模式，避免大量内存复制
void* setup_direct_access_mode(void *shm_addr) {
    debug_log("Setting up direct access mode for shared memory at %p", shm_addr);
    
    // 获取头部信息
    shm_header_t *header = (shm_header_t *)shm_addr;
    
    // 验证魔数
    if (header->magic_number != SHM_MAGIC_NUMBER) {
        debug_log("ERROR: Invalid magic number in shared memory: %x (expected %x)", 
                  header->magic_number, SHM_MAGIC_NUMBER);
        return NULL;
    }
    
    // 获取共享内存中的各部分地址
    char *readnames = (char *)shm_addr + header->readnames_offset;
    bam_read_idx_record *records = (bam_read_idx_record *)((char *)shm_addr + header->records_offset);
    bam_read_idx *bri = (bam_read_idx *)((char *)shm_addr + header->bri_offset);
    
    debug_log("Shared memory layout: bri=%p, readnames=%p, records=%p", 
              bri, readnames, records);
    debug_log("Header info: name_count=%zu, record_count=%zu", 
              header->name_count_bytes, header->record_count);
    
    // 创建上下文信息
    client_context_t *ctx = (client_context_t *)malloc(sizeof(client_context_t));
    if (!ctx) {
        debug_log("ERROR: Failed to allocate memory for client context");
        return NULL;
    }
    
    // 初始化上下文
    ctx->shm_addr = shm_addr;
    ctx->record_count = header->record_count;
    ctx->readnames = readnames;
    ctx->records = records;
    
    // 设置bri结构
    bri->readnames = readnames;
    bri->name_count_bytes = header->name_count_bytes;
    bri->record_count = header->record_count;
    
    // 重要：将bri->records指向我们的上下文，这样search函数可以通过它来访问
    bri->records = (bam_read_idx_record*)ctx;
    
    // 修改header标记使用直接访问模式
    header->direct_access_mode = 1;
    
    debug_log("Direct access mode setup completed");
    return ctx;
}

// 修改连接到共享内存的函数，使用直接访问模式
static PyObject* py_connect_to_server(PyObject *self, PyObject *args)
{
    const char *bam_file;
    
    if (!PyArg_ParseTuple(args, "s", &bam_file)) {
        return NULL;
    }
    
    debug_log("Connecting to server for BAM file: %s", bam_file);
    
    // 检查密钥文件和大小文件是否存在
    if (access(SHM_KEY_PATH, F_OK) == -1 || access(SHM_SIZE_PATH, F_OK) == -1) {
        debug_log("Key file or size file does not exist, server not running");
        Py_RETURN_NONE;
    }
    
    // 读取共享内存大小
    FILE *size_fp = fopen(SHM_SIZE_PATH, "r");
    if (size_fp == NULL) {
        debug_log("Failed to open size file: %s", strerror(errno));
        Py_RETURN_NONE;
    }
    
    size_t total_size;
    if (fscanf(size_fp, "%zu", &total_size) != 1) {
        debug_log("Failed to read size from file");
        fclose(size_fp);
        Py_RETURN_NONE;
    }
    fclose(size_fp);
    
    debug_log("Read shared memory size: %zu bytes", total_size);
    
    // 生成共享内存键
    key_t key = ftok(SHM_KEY_PATH, SHM_PROJ_ID);
    if (key == -1) {
        debug_log("Failed to generate key with ftok: %s", strerror(errno));
        Py_RETURN_NONE;
    }
    
    debug_log("Generated shared memory key: %d", key);
    
    // 获取共享内存段
    int shm_id = shmget(key, total_size, 0666);
    if (shm_id == -1) {
        debug_log("Failed to get shared memory: %s", strerror(errno));
        Py_RETURN_NONE;
    }
    
    debug_log("Found shared memory segment: ID=%d", shm_id);
    
    // 附加到共享内存段
    void *shm_addr = shmat(shm_id, NULL, 0);
    if (shm_addr == (void *)-1) {
        debug_log("Failed to attach to shared memory: %s", strerror(errno));
        Py_RETURN_NONE;
    }
    
    debug_log("Attached to shared memory at %p", shm_addr);
    
    // 获取头部信息
    shm_header_t *header = (shm_header_t *)shm_addr;
    
    // 验证魔数
    if (header->magic_number != SHM_MAGIC_NUMBER) {
        debug_log("Invalid magic number: %x (expected %x)", header->magic_number, SHM_MAGIC_NUMBER);
        shmdt(shm_addr);
        Py_RETURN_NONE;
    }
    
    debug_log("Magic number verified");
    
    // 验证BAM文件路径是否匹配
    if (strcmp(header->bam_path, bam_file) != 0) {
        debug_log("BAM file mismatch: '%s' vs '%s'", header->bam_path, bam_file);
        shmdt(shm_addr);
        Py_RETURN_NONE;
    }
    
    debug_log("BAM file path verified");
    
    // 在验证共享内存结构之前，使用直接访问模式
    void *client_ctx = setup_direct_access_mode(shm_addr);
    
    if (!client_ctx) {
        debug_log("ERROR: Failed to setup direct access mode");
        shmdt(shm_addr);
        Py_RETURN_NONE;
    }
    
    // 获取bri结构体指针
    bam_read_idx *bri = (bam_read_idx *)((char *)shm_addr + header->bri_offset);
    
    // 验证设置是否成功
    debug_log("Validating direct access setup...");
    
    // 选择几个示例记录点进行验证
    size_t check_points[] = {0, 1, 2, 3, 4, 100000, 1000000};
    int validation_passed = 1;
    
    for (size_t i = 0; i < sizeof(check_points)/sizeof(check_points[0]); i++) {
        size_t idx = check_points[i];
        if (idx >= header->record_count) continue;
        
        client_context_t *ctx = (client_context_t*)bri->records;
        const char *name = get_record_readname(ctx, idx);
        
        if (!name) {
            debug_log("ERROR: Failed to get readname for record %zu", idx);
            validation_passed = 0;
            break;
        }
        
        debug_log("Record %zu validation: name=%s", idx, name);
    }
    
    if (!validation_passed) {
        debug_log("CRITICAL: Direct access validation failed, detaching from shared memory");
        free(client_ctx);
        shmdt(shm_addr);
        Py_RETURN_NONE;
    } else {
        debug_log("Direct access validation passed!");
    }
    
    // 创建返回值
    PyObject *shm_addr_obj = PyLong_FromVoidPtr(shm_addr);
    PyObject *bri_ptr_obj = PyLong_FromVoidPtr(bri);
    PyObject *ctx_obj = PyLong_FromVoidPtr(client_ctx);
    
    if (!shm_addr_obj || !bri_ptr_obj || !ctx_obj) {
        debug_log("ERROR: Failed to create Python objects");
        free(client_ctx);
        shmdt(shm_addr);
        Py_XDECREF(shm_addr_obj);
        Py_XDECREF(bri_ptr_obj);
        Py_XDECREF(ctx_obj);
        Py_RETURN_NONE;
    }
    
    // 返回三个值：共享内存地址、bri指针和客户端上下文
    PyObject *result = PyTuple_New(3);
    PyTuple_SetItem(result, 0, shm_addr_obj);
    PyTuple_SetItem(result, 1, bri_ptr_obj);
    PyTuple_SetItem(result, 2, ctx_obj);
    
    debug_log("Connection successful with direct access mode");
    
    return result;
}

// 断开与共享内存的连接
static PyObject* py_disconnect_from_server(PyObject *self, PyObject *args)
{
    PyObject *shm_addr_obj;
    PyObject *ctx_obj = NULL;
    
    if (!PyArg_ParseTuple(args, "O|O", &shm_addr_obj, &ctx_obj)) {
        return NULL;
    }
    
    void *shm_addr = PyLong_AsVoidPtr(shm_addr_obj);
    if (shm_addr == NULL && PyErr_Occurred()) {
        return NULL;
    }
    
    debug_log("Disconnecting from server, shm_addr=%p", shm_addr);
    
    // 获取BRI结构
    shm_header_t *header = (shm_header_t *)shm_addr;
    bam_read_idx *bri = (bam_read_idx *)((char *)shm_addr + header->bri_offset);
    
    // 释放客户端上下文
    if (ctx_obj != NULL) {
        void *client_ctx = PyLong_AsVoidPtr(ctx_obj);
        if (client_ctx != NULL) {
            debug_log("Freeing client context at %p", client_ctx);
            
            // 检查bri->records是否指向客户端上下文（在direct_access模式下是这种情况）
            if (bri && bri->records == (bam_read_idx_record*)client_ctx) {
                debug_log("Detected that bri->records points to client context, setting to NULL to avoid double free");
                bri->records = NULL; // 避免double free
            }
            
            free(client_ctx);
        }
    }
    
    // 释放本地记录数组，但只有在它不是客户端上下文的情况下
    if (bri && bri->records) {
        // 检查bri->records是否指向共享内存外部的地址
        bam_read_idx_record *shm_records = (bam_read_idx_record *)((char *)shm_addr + header->records_offset);
        if (bri->records != shm_records) {
            debug_log("Freeing local records array at %p", bri->records);
            free(bri->records);
            bri->records = NULL;
        }
    }
    
    // 从共享内存分离
    if (shmdt(shm_addr) == -1) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to detach from shared memory");
        debug_log("Failed to detach from shared memory: %s", strerror(errno));
        return NULL;
    }
    
    debug_log("Successfully disconnected from server");
    Py_RETURN_NONE;
}

// 验证共享内存中的结构体状态
void validate_shm_structures(void *shm_addr) {
    shm_header_t *header = (shm_header_t *)shm_addr;
    bam_read_idx *bri = (bam_read_idx *)((char *)shm_addr + header->bri_offset);
    
    debug_log("Validating shared memory structures");
    debug_log("Header stats: magic=%x, name_count=%zu, record_count=%zu", 
              header->magic_number, header->name_count_bytes, header->record_count);
    debug_log("BRI stats: name_count=%zu, record_count=%zu, readnames=%p, records=%p", 
              bri->name_count_bytes, bri->record_count, bri->readnames, bri->records);
    
    // 检查前5条记录
    bam_read_idx_record *records = bri->records;
    for (size_t i = 0; i < 5 && i < bri->record_count; i++) {
        debug_log("Record %zu: name_ptr=%p, file_offset=%zu", 
                 i, records[i].read_name.ptr, records[i].file_offset);
        if (records[i].read_name.ptr != NULL) {
            debug_log("  Name content: '%s'", records[i].read_name.ptr);
        } else {
            debug_log("  WARNING: NULL name pointer");
        }
    }
}

// 创建共享内存并复制BRI索引到共享内存中
static PyObject* py_create_server(PyObject *self, PyObject *args)
{
    const char *bam_file;
    const char *index_file;
    
    if (!PyArg_ParseTuple(args, "ss", &bam_file, &index_file)) {
        return NULL;
    }

    debug_log("Creating server for BAM: %s, Index: %s", bam_file, index_file);
    
    // 加载BRI索引
    bam_read_idx *bri = initialize_bri(bam_file, index_file);
    if (!bri) {
        PyErr_SetString(PyExc_RuntimeError, "Could not initialize bri.");
        return NULL;
    }
    
    debug_log("BRI loaded with %zu records, %zu name bytes", bri->record_count, bri->name_count_bytes);
    
    // 计算需要的共享内存总大小
    size_t total_size = calculate_total_size(bri);
    
    // 创建元数据文件，存储共享内存大小
    FILE *size_fp = fopen(SHM_SIZE_PATH, "w");
    if (size_fp == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create size file");
        bam_read_idx_destroy(bri);
        return NULL;
    }
    fprintf(size_fp, "%zu", total_size);
    fclose(size_fp);
    
    // 创建密钥文件
    FILE *key_fp = fopen(SHM_KEY_PATH, "w");
    if (key_fp == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create key file");
        bam_read_idx_destroy(bri);
        unlink(SHM_SIZE_PATH);
        return NULL;
    }
    fclose(key_fp);
    
    // 生成共享内存键
    key_t key = ftok(SHM_KEY_PATH, SHM_PROJ_ID);
    if (key == -1) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to generate key with ftok");
        debug_log("ftok error: %s", strerror(errno));
        bam_read_idx_destroy(bri);
        unlink(SHM_KEY_PATH);
        unlink(SHM_SIZE_PATH);
        return NULL;
    }
    
    debug_log("Generated shared memory key: %d", key);
    
    // 创建共享内存段
    int shm_id = shmget(key, total_size, IPC_CREAT | 0666);
    if (shm_id == -1) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create shared memory");
        debug_log("shmget error: %s", strerror(errno));
        bam_read_idx_destroy(bri);
        unlink(SHM_KEY_PATH);
        unlink(SHM_SIZE_PATH);
        return NULL;
    }
    
    debug_log("Created shared memory segment: ID=%d, size=%zu", shm_id, total_size);
    
    // 附加到共享内存段
    void *shm_addr = shmat(shm_id, NULL, 0);
    if (shm_addr == (void *)-1) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to attach to shared memory");
        debug_log("shmat error: %s", strerror(errno));
        bam_read_idx_destroy(bri);
        shmctl(shm_id, IPC_RMID, NULL);
        unlink(SHM_KEY_PATH);
        unlink(SHM_SIZE_PATH);
        return NULL;
    }
    
    debug_log("Attached to shared memory at %p", shm_addr);
    
    // 初始化整个共享内存为0，确保没有垃圾数据
    memset(shm_addr, 0, total_size);
    
    // 初始化共享内存头部
    shm_header_t *header = (shm_header_t *)shm_addr;
    header->total_size = total_size;
    header->bri_offset = sizeof(shm_header_t);
    header->readnames_offset = header->bri_offset + sizeof(bam_read_idx);
    header->records_offset = header->readnames_offset + bri->name_count_bytes;
    header->name_count_bytes = bri->name_count_bytes;
    header->record_count = bri->record_count;
    header->magic_number = SHM_MAGIC_NUMBER;  // 设置魔数
    header->new_approach = 1;  // 标记使用新方法
    strncpy(header->bam_path, bam_file, sizeof(header->bam_path) - 1);
    header->bam_path[sizeof(header->bam_path) - 1] = '\0';
    
    debug_log("Initialized header: bri_offset=%zu, readnames_offset=%zu, records_offset=%zu", 
              header->bri_offset, header->readnames_offset, header->records_offset);
    
    // 复制bri结构体到共享内存
    bam_read_idx *shm_bri = (bam_read_idx *)((char *)shm_addr + header->bri_offset);
    memcpy(shm_bri, bri, sizeof(bam_read_idx));
    debug_log("Copied BRI structure to shared memory at %p", shm_bri);
    
    // 复制readnames数据到共享内存
    char *shm_readnames = (char *)shm_addr + header->readnames_offset;
    memcpy(shm_readnames, bri->readnames, bri->name_count_bytes);
    debug_log("Copied readnames data to shared memory at %p, bytes=%zu", 
              shm_readnames, bri->name_count_bytes);
    
    // 检查前几个字符串是否正确复制
    debug_log("First few readnames: '%s', '%s'", 
              shm_readnames, shm_readnames + strlen(shm_readnames) + 1);
    
    // 复制records数据到共享内存，但先保留原始offset
    bam_read_idx_record *shm_records = (bam_read_idx_record *)((char *)shm_addr + header->records_offset);
    memcpy(shm_records, bri->records, bri->record_count * sizeof(bam_read_idx_record));
    debug_log("Copied %zu records to shared memory at %p", bri->record_count, shm_records);
    
    // 把records中的read_name.ptr转换为offset存储
    for (size_t i = 0; i < bri->record_count; i++) {
        // 计算在原始内存中的偏移量
        const char *ptr = bri->records[i].read_name.ptr;
        if (ptr != NULL) {
            // 计算字符偏移量
            size_t offset = (size_t)(ptr - bri->readnames);
            
            // 验证偏移量在合理范围内
            if (offset < bri->name_count_bytes) {
                // 彻底解决方案：不再使用union，直接在共享内存中只存储偏移量
                // 完全清零整个记录结构
                memset(&shm_records[i], 0, sizeof(bam_read_idx_record));
                
                // 直接将偏移量值复制到read_name的内存位置
                memcpy(&shm_records[i].read_name, &offset, sizeof(size_t));
                
                // 日志记录验证偏移量是否正确 (只记录少量样本)
                if (i < 5 || i == bri->record_count - 1 || i % 10000000 == 0) {
                    debug_log("Record %zu: Original ptr=%p, offset=%zu, name=%s", 
                             i, ptr, offset, ptr);
                    debug_log("  Verification: Name at offset in shm: %s", 
                             shm_readnames + offset);
                    
                    // 验证内存
                    size_t stored_offset;
                    memcpy(&stored_offset, &shm_records[i].read_name, sizeof(size_t));
                    debug_log("  Stored offset value (direct memory): %zu", stored_offset);
                    
                    // 输出详细的二进制内容
                    unsigned char *bytes = (unsigned char*)&shm_records[i].read_name;
                    debug_log("  Memory content: %02x %02x %02x %02x %02x %02x %02x %02x",
                              bytes[0], bytes[1], bytes[2], bytes[3], 
                              bytes[4], bytes[5], bytes[6], bytes[7]);
                }
                
                // 复制文件偏移量
                shm_records[i].file_offset = bri->records[i].file_offset;
            } else {
                debug_log("ERROR: Invalid offset calculation for record %zu: %zu (exceeds %zu), ptr=%p, readnames=%p", 
                         i, offset, bri->name_count_bytes, ptr, bri->readnames);
                // 安全措施
                memset(&shm_records[i], 0, sizeof(bam_read_idx_record));
                size_t zero = 0;
                memcpy(&shm_records[i].read_name, &zero, sizeof(size_t));
                shm_records[i].file_offset = bri->records[i].file_offset;
            }
        } else {
            debug_log("ERROR: NULL name pointer for record %zu", i);
            // 安全措施
            memset(&shm_records[i], 0, sizeof(bam_read_idx_record));
            size_t zero = 0;
            memcpy(&shm_records[i].read_name, &zero, sizeof(size_t));
            shm_records[i].file_offset = bri->records[i].file_offset;
        }
    }
    
    // 重建共享内存中的指针关系
    // rebuild_pointers_in_shm(shm_addr);
    
    // 验证共享内存中的结构
    validate_shm_structures(shm_addr);
    
    // 从共享内存分离
    if (shmdt(shm_addr) == -1) {
        debug_log("WARNING: Failed to detach from shared memory: %s", strerror(errno));
    }
    
    // 释放原始bri结构体
    bam_read_idx_destroy(bri);
    debug_log("Server setup completed, shared memory ID: %d", shm_id);
    
    return Py_BuildValue("i", shm_id);
}

// 关闭服务器并释放资源
static PyObject* py_stop_server(PyObject *self, PyObject *args)
{
    debug_log("Stopping server");
    
    // 如果文件不存在，不要立即报错，尝试清理
    int missing_files = 0;
    
    if (access(SHM_KEY_PATH, F_OK) == -1 || access(SHM_SIZE_PATH, F_OK) == -1) {
        missing_files = 1;
        debug_log("Key file or size file does not exist");
    }
    
    size_t total_size = 0;
    key_t key = 0;
    int shm_id = -1;
    
    // 如果文件存在，获取共享内存信息
    if (!missing_files) {
        // 读取共享内存大小
        FILE *size_fp = fopen(SHM_SIZE_PATH, "r");
        if (size_fp == NULL) {
            missing_files = 1;
            debug_log("Failed to open size file: %s", strerror(errno));
        } else {
            if (fscanf(size_fp, "%zu", &total_size) != 1) {
                fclose(size_fp);
                missing_files = 1;
                debug_log("Failed to read size from file");
            } else {
                fclose(size_fp);
                debug_log("Read shared memory size: %zu bytes");
            }
        }
        
        // 获取共享内存键和ID
        if (!missing_files) {
            key = ftok(SHM_KEY_PATH, SHM_PROJ_ID);
            if (key == -1) {
                missing_files = 1;
                debug_log("Failed to generate key with ftok: %s", strerror(errno));
            } else {
                shm_id = shmget(key, total_size, 0666);
                if (shm_id == -1) {
                    missing_files = 1;
                    debug_log("Failed to get shared memory: %s", strerror(errno));
                } else {
                    debug_log("Found shared memory segment: ID=%d", shm_id);
                }
            }
        }
    }
    
    // 尝试删除共享内存段
    if (!missing_files && shm_id != -1) {
        if (shmctl(shm_id, IPC_RMID, NULL) == -1) {
            // 如果删除失败，记录但继续运行
            PyErr_WarnFormat(PyExc_RuntimeWarning, 1, 
                            "Failed to remove shared memory segment: %s", strerror(errno));
            debug_log("Failed to remove shared memory segment: %s", strerror(errno));
        } else {
            debug_log("Removed shared memory segment ID=%d", shm_id);
        }
    }
    
    // 删除密钥文件和大小文件
    int key_removed = 0;
    int size_removed = 0;
    
    if (access(SHM_KEY_PATH, F_OK) != -1) {
        if (unlink(SHM_KEY_PATH) != -1) {
            key_removed = 1;
            debug_log("Removed key file: %s", SHM_KEY_PATH);
        } else {
            debug_log("Failed to remove key file: %s", strerror(errno));
        }
    }
    
    if (access(SHM_SIZE_PATH, F_OK) != -1) {
        if (unlink(SHM_SIZE_PATH) != -1) {
            size_removed = 1;
            debug_log("Removed size file: %s", SHM_SIZE_PATH);
        } else {
            debug_log("Failed to remove size file: %s", strerror(errno));
        }
    }
    
    // 如果我们没有成功清理任何东西，报告错误
    if (missing_files && !key_removed && !size_removed) {
        PyErr_SetString(PyExc_RuntimeError, "Server not running or already stopped");
        debug_log("Server not running or already stopped");
        return NULL;
    }
    
    debug_log("Server stopped successfully");
    Py_RETURN_NONE;
}

static PyObject* py_query_by_readname(PyObject *self, PyObject *args)
{
    unsigned long long bri_ptr;
    PyObject *read_names_list;
    const char *input_bam;
    if (!PyArg_ParseTuple(args, "sKO", &input_bam, &bri_ptr, &read_names_list)) {
        return NULL;
    }
    
    if (!PyList_Check(read_names_list)) {
        PyErr_SetString(PyExc_TypeError, "read_names must be a list.");
        return NULL;
    }

    bam_read_idx *bri = (bam_read_idx *) bri_ptr;
    
    debug_log("Query by readname using BRI at %p for BAM %s", bri, input_bam);
    debug_log("BRI stats: readnames=%p, records=%p, name_count=%zu, record_count=%zu", 
              bri->readnames, bri->records, bri->name_count_bytes, bri->record_count);
    
    // 检查BRI指针的有效性
    if (bri == NULL || bri->readnames == NULL || bri->records == NULL || 
        bri->record_count == 0 || bri->name_count_bytes == 0) {
        debug_log("ERROR: Invalid BRI structure");
        PyErr_SetString(PyExc_RuntimeError, "Invalid BRI structure");
        return NULL;
    }
    
    Py_ssize_t size = PyList_Size(read_names_list);
    char **read_names = malloc(size * sizeof(char *));
    
    debug_log("Processing %zd read names", size);
    
    for (Py_ssize_t i = 0; i < size; i++) {
        PyObject *item = PyList_GetItem(read_names_list, i);
        if (!PyUnicode_Check(item)) {
            free(read_names);
            PyErr_SetString(PyExc_TypeError, "All read names must be strings.");
            return NULL;
        }
        
        read_names[i] = (char*)PyUnicode_AsUTF8(item);
        if (i < 5) {  // 只记录前5个，避免日志过大
            debug_log("Read name %zd: %s", i, read_names[i]);
        }
    }

    kstring_t ks = {0, 0, NULL};  // 初始化 kstring_t

    debug_log("Calling query_by_readname function");
    query_by_readname(input_bam, bri, read_names, (uint64_t)size, &ks);
    debug_log("Query completed, result size: %zu bytes", ks.l);
    
    if (ks.s == NULL) {
        debug_log("ERROR: NULL result from query_by_readname");
        free(read_names);
        PyErr_SetString(PyExc_RuntimeError, "Failed to get results from BAM file");
        return NULL;
    }

    PyObject *result = Py_BuildValue("s", ks.s);  // 将 kstring_t 转换为 Python 字符串
    if (result == NULL) {
        debug_log("ERROR: Failed to convert result to Python string");
    }

    free(ks.s);  // 释放 kstring_t
    free(read_names);  // 释放 read_names

    debug_log("Function completed successfully");
    return result;  // 返回 Python 字符串
}

static PyMethodDef BriMethods[] = {
    {"initialize_bri",  py_initialize_bri, METH_VARARGS, "Initialize bam_read_idx."},
    {"query_by_readname",  py_query_by_readname, METH_VARARGS, "Query by read name."},
    {"create_server",  py_create_server, METH_VARARGS, "Start server mode and share BRI index."},
    {"connect_to_server",  py_connect_to_server, METH_VARARGS, "Connect to running server."},
    {"disconnect_from_server",  py_disconnect_from_server, METH_VARARGS, "Disconnect from server."},
    {"stop_server",  py_stop_server, METH_VARARGS, "Stop server and free resources."},
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
