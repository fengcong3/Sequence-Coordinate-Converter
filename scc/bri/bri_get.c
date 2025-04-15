//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// bri - simple utility to provide random access to
//       bam records by read name
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "bri_get.h"
//
// Getopt
//
enum {
    OPT_HELP = 1,
};

static const char* shortopts = ":i:"; // placeholder
static const struct option longopts[] = {
    { "help",                      no_argument,       NULL, OPT_HELP },
    { "index",               required_argument,       NULL,      'i' },
    { NULL, 0, NULL, 0 }
};

void print_usage_get()
{
    fprintf(stderr, "usage: bri get [-i <index_filename.bri>] <input.bam> <readname>\n");
}

// comparator used by bsearch, direct strcmp through the name pointer
int compare_records_by_readname_ptr(const void* r1, const void* r2)
{
    const char* n1 = ((bam_read_idx_record*)r1)->read_name.ptr;
    const char* n2 = ((bam_read_idx_record*)r2)->read_name.ptr;
    return strcmp(n1, n2);
}

//
void bam_read_idx_get_range(const bam_read_idx* bri, const char* readname, bam_read_idx_record** start, bam_read_idx_record** end)
{
    // construct a query record to pass to bsearch
    bam_read_idx_record query;
    query.read_name.ptr = readname;
    query.file_offset = 0;

    // if rec is NULL then readname does not appear in index
    bam_read_idx_record* rec = 
        bsearch(&query, bri->records, bri->record_count, sizeof(bam_read_idx_record), compare_records_by_readname_ptr);
    if(rec == NULL) {
        *start = NULL;
        *end = NULL;
        return;
    }

    // rec points to a valid record, but it can be an arbitrary record in the range
    // move start to the first record in the range, and end to be one past the end
    size_t sri, eri;
    sri = eri = rec - bri->records;
    assert(bri->records[sri].file_offset == rec->file_offset);

    while(sri > 0 && bri->records[sri].read_name.ptr == bri->records[sri - 1].read_name.ptr) {
        sri -= 1;
    }
    assert(strcmp(bri->records[sri].read_name.ptr, readname) == 0);

    do {
        eri += 1;
    } while(eri < bri->record_count && bri->records[eri].read_name.ptr == bri->records[sri].read_name.ptr);
    assert(eri == bri->record_count || strcmp(bri->records[eri].read_name.ptr, readname) != 0);
    //fprintf(stderr, "r: %zu sri: %zu eri: %zu\n", rec - bri->records, sri, eri);
    *start = &bri->records[sri];
    *end = &bri->records[eri];
}

//
void bam_read_idx_get_by_record(htsFile* fp, bam_hdr_t* hdr, bam1_t* b, bam_read_idx_record* bri_record)
{
    int ret = bgzf_seek(fp->fp.bgzf, bri_record->file_offset, SEEK_SET);
    if(ret != 0) {
        fprintf(stderr, "[bri] bgzf_seek failed\n");
        exit(EXIT_FAILURE);
    }

    ret = sam_read1(fp, hdr, b);
    if(ret < 0) {
        fprintf(stderr, "[bri] sam_read1 failed\n");
        exit(EXIT_FAILURE);
    }
}

// Initialize bam_read_idx
bam_read_idx* initialize_bri(const char* bam_file, const char* index_file) {
    return bam_read_idx_load(bam_file, index_file);
}

// Query by readname --fengcong
// void query_by_readname(const char* input_bam, bam_read_idx* bri,
//     char** read_names, uint64_t read_name_count, const char* output_sam){
//     htsFile* bam_fp = hts_open(input_bam, "r");
//     bam_hdr_t* h = sam_hdr_read(bam_fp);
//     htsFile* out_fp = hts_open(output_sam, "w");

//     bam_read_idx_record* start;
//     bam_read_idx_record* end;

//     //write header
    
//     int retx = sam_hdr_write(out_fp, h);//写入到out_fp中,header 和 b
//     if(retx < 0) {
//         fprintf(stderr, "[bri] sam_hdr_write failed\n");
//         exit(EXIT_FAILURE);
//     }

//     for(int i = 0; i < read_name_count; i++) {
//         char* readname = read_names[i]; // 一个readname
//         bam_read_idx_get_range(bri, readname, &start, &end); // 根据一个readname 获得一个start end
//         bam1_t* b = bam_init1();
//         while(start != end) { // 同样的readname 可能有两条记录
//             bam_read_idx_get_by_record(bam_fp, h, b, start); //得到的record写入到b中
//             int ret = sam_write1(out_fp, h, b);//写入到out_fp中,header 和 b
//             if(ret < 0) {
//                 fprintf(stderr, "[bri] sam_write1 failed\n");
//                 exit(EXIT_FAILURE);
//             }

//             start++;
//         }
//         bam_destroy1(b);
//     }

//     hts_close(out_fp);
//     bam_hdr_destroy(h);
//     hts_close(bam_fp);
//     // bam_read_idx_destroy(bri);
//     // bri = NULL;
// }

void query_by_readname(const char* input_bam, bam_read_idx* bri,
    char** read_names, uint64_t read_name_count, kstring_t *ks_out)
{
    // 添加安全检查
    if (bri == NULL) {
        fprintf(stderr, "[bri] Error: NULL bri pointer passed to query_by_readname\n");
        return;
    }
    
    // 检查是否在直接访问模式
    int direct_mode = 0;
    client_context_t *ctx = NULL;
    
    // 如果records指向了一个client_context_t结构体，那么我们是在直接访问模式下
    if (bri->records != NULL && bri->record_count > 0) {
        // 尝试获取readname以测试是否在直接访问模式
        ctx = (client_context_t*)bri->records;
        if (get_record_readname(ctx, 0) != NULL) {
            direct_mode = 1;
        }
    }
    
    if (!direct_mode && (bri->readnames == NULL || bri->record_count == 0)) {
        fprintf(stderr, "[bri] Error: Invalid bri structure (NULL readnames or record_count = 0)\n");
        return;
    }
    
    htsFile* bam_fp = hts_open(input_bam, "r");
    if (bam_fp == NULL) {
        fprintf(stderr, "[bri] Error: Could not open BAM file: %s\n", input_bam);
        return;
    }
    
    bam_hdr_t* h = sam_hdr_read(bam_fp);
    if (h == NULL) {
        fprintf(stderr, "[bri] Error: Could not read header from BAM file\n");
        hts_close(bam_fp);
        return;
    }

    // 确保 ks_out 已初始化
    ks_clear(ks_out);
    // 将整个 header 写入 kstring
    kputs(h->text, ks_out);

    for(uint64_t i = 0; i < read_name_count; i++) {
        char* readname = read_names[i];
        
        // 添加安全检查
        if (readname == NULL) {
            fprintf(stderr, "[bri] Warning: NULL readname at index %lu, skipping\n", i);
            continue;
        }
        
        // 根据模式选择不同的查找方式
        if (direct_mode) {
            size_t start_idx, end_idx;
            if (!direct_search_record_by_name(ctx, readname, &start_idx, &end_idx)) {
                // 找不到此readname，跳过
                continue;
            }
            
            bam1_t* b = bam_init1();
            if (b == NULL) {
                fprintf(stderr, "[bri] Error: Failed to allocate memory for bam1_t\n");
                break;
            }
            
            for (size_t idx = start_idx; idx < end_idx; idx++) {
                // 获取文件偏移量
                size_t file_offset = ctx->records[idx].file_offset;
                
                // 定位并读取BAM记录
                int ret = bgzf_seek(bam_fp->fp.bgzf, file_offset, SEEK_SET);
                if (ret != 0) {
                    fprintf(stderr, "[bri] bgzf_seek failed\n");
                    continue;
                }
                
                ret = sam_read1(bam_fp, h, b);
                if (ret < 0) {
                    fprintf(stderr, "[bri] sam_read1 failed\n");
                    continue;
                }
                
                kstring_t str = {0, 0, NULL};
                ret = sam_format1(h, b, &str);
                if (ret < 0 || str.s == NULL) {
                    fprintf(stderr, "[bri] Error: Failed to format record\n");
                    free(str.s);
                    continue;
                }
                
                kputs(str.s, ks_out);
                kputc('\n', ks_out);
                free(str.s);
            }
            bam_destroy1(b);
        } else {
            // 原始的获取记录范围的方法
            bam_read_idx_record* start;
            bam_read_idx_record* end;
            bam_read_idx_get_range(bri, readname, &start, &end);
            
            if (start == NULL || end == NULL) {
                // 找不到此readname，跳过
                continue;
            }
            
            bam1_t* b = bam_init1();
            if (b == NULL) {
                fprintf(stderr, "[bri] Error: Failed to allocate memory for bam1_t\n");
                break;
            }
            
            while(start != end) {
                // 检查start是否为有效指针
                if (start < bri->records || start >= bri->records + bri->record_count) {
                    fprintf(stderr, "[bri] Error: Invalid record pointer\n");
                    break;
                }
                
                bam_read_idx_get_by_record(bam_fp, h, b, start);
                
                kstring_t str = {0, 0, NULL};
                int ret = sam_format1(h, b, &str);
                if (ret < 0 || str.s == NULL) {
                    fprintf(stderr, "[bri] Error: Failed to format record\n");
                    free(str.s);
                    continue;
                }
                
                kputs(str.s, ks_out);
                kputc('\n', ks_out);
                free(str.s);

                start++;
            }
            bam_destroy1(b);
        }
    }

    bam_hdr_destroy(h);
    hts_close(bam_fp);
}

//
int bam_read_idx_get_main(int argc, char** argv)
{
    char* input_bri = NULL; // bri file name/path

    int die = 0;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        switch (c) {
            case OPT_HELP:
                print_usage_get();
                exit(EXIT_SUCCESS);
            case 'i':
                input_bri = optarg;
        }
    }
    
    if (argc - optind < 2) {
        fprintf(stderr, "bri get: not enough arguments\n");
        die = 1;
    }

    if(die) {
        print_usage_get();
        exit(EXIT_FAILURE);
    }

    char* input_bam = argv[optind++];

    bam_read_idx* bri = bam_read_idx_load(input_bam, input_bri);
    
    htsFile* bam_fp = hts_open(input_bam, "r");
    bam_hdr_t* h = sam_hdr_read(bam_fp);
    htsFile* out_fp = hts_open("-", "w");

    bam_read_idx_record* start;
    bam_read_idx_record* end;

    for(int i = optind; i < argc; i++) {
        char* readname = argv[i]; // 一个readname
        bam_read_idx_get_range(bri, readname, &start, &end); // 根据一个readname 获得一个start end
        bam1_t* b = bam_init1();
        while(start != end) { // 同样的readname 可能有两条记录
            
            bam_read_idx_get_by_record(bam_fp, h, b, start); //得到的record写入到b中
            int ret = sam_write1(out_fp, h, b);//写入到out_fp中,header 和 b
            if(ret < 0) {
                fprintf(stderr, "[bri] sam_write1 failed\n");
                exit(EXIT_FAILURE);
            }

            start++;
        }
        bam_destroy1(b);
    }

    hts_close(out_fp);
    bam_hdr_destroy(h);
    hts_close(bam_fp);
    bam_read_idx_destroy(bri);
    bri = NULL;

    return 0;
}

// 从客户端上下文中获取指定索引的readname
const char* get_record_readname(client_context_t *ctx, size_t index) {
    if (index >= ctx->record_count) {
        return NULL;
    }
    
    // 从共享内存读取偏移量
    size_t offset;
    memcpy(&offset, &ctx->records[index].read_name, sizeof(size_t));
    
    // 验证偏移量
    if (offset >= 0) {
        return ctx->readnames + offset;
    }
    
    return NULL;
}

// 二分搜索找到记录，使用按需访问模式
int direct_search_record_by_name(client_context_t *ctx, const char *readname, size_t *start_idx, size_t *end_idx) {
    // 如果没有记录，直接返回失败
    if (ctx->record_count == 0) {
        return 0;
    }
    
    // 二分查找
    size_t left = 0;
    size_t right = ctx->record_count - 1;
    size_t mid;
    int found = 0;
    
    while (left <= right) {
        mid = left + (right - left) / 2;
        const char *current_name = get_record_readname(ctx, mid);
        
        if (!current_name) {
            return 0; // 出现错误
        }
        
        int cmp = strcmp(readname, current_name);
        if (cmp == 0) {
            found = 1;
            break;
        } else if (cmp < 0) {
            if (mid == 0) break;
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }
    
    if (!found) {
        return 0;
    }
    
    // 找到匹配项后，向前和向后扩展查找范围
    size_t first = mid;
    while (first > 0) {
        const char *prev_name = get_record_readname(ctx, first - 1);
        if (!prev_name || strcmp(readname, prev_name) != 0) {
            break;
        }
        first--;
    }
    
    size_t last = mid;
    while (last < ctx->record_count - 1) {
        const char *next_name = get_record_readname(ctx, last + 1);
        if (!next_name || strcmp(readname, next_name) != 0) {
            break;
        }
        last++;
    }
    
    *start_idx = first;
    *end_idx = last + 1; // 返回最后一个匹配项的下一个位置
    return 1;
}
