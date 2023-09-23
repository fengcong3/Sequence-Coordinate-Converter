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
    char** read_names, uint64_t read_name_count, kstring_t *ks_out){
    htsFile* bam_fp = hts_open(input_bam, "r");
    bam_hdr_t* h = sam_hdr_read(bam_fp);

    bam_read_idx_record* start;
    bam_read_idx_record* end;

    // 确保 ks_out 已初始化
    ks_clear(ks_out);
    //write header
    
    // 将整个 header 写入 kstring
    kputs(h->text, ks_out);

    for(int i = 0; i < read_name_count; i++) {
        char* readname = read_names[i]; // 一个readname
        bam_read_idx_get_range(bri, readname, &start, &end); // 根据一个readname 获得一个start end
        bam1_t* b = bam_init1();
        while(start != end) { // 同样的readname 可能有两条记录
            bam_read_idx_get_by_record(bam_fp, h, b, start); //得到的record写入到b中
            kstring_t str = {0, 0, NULL};
            sam_format1(h, b, &str);
            kputs(str.s, ks_out);
            kputc('\n', ks_out);
            free(str.s);

            start++;
        }
        bam_destroy1(b);
    }

    // hts_close(out_fp);
    bam_hdr_destroy(h);
    hts_close(bam_fp);
    // bam_read_idx_destroy(bri);
    // bri = NULL;
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
