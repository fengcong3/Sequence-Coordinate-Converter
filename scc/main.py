#coding:utf-8
# Author: fengcong@caas.cn
# Date: 2023-09-22
# SequenceCoordinateConverter:
# --------|---------reference 1----------------- chr1A
#     ----|-read1----->
#     <---|--read2-----
#       --|---read3----->

## "|" represents the pos need to be converted

# -----------------reference 2------|----------- chr1A
#                               ----|-read1----->
#                               <---|--read2-----
#                                 --|---read3----->

## "|" represents the pos has been converted

### So, When we have a position (Pos1) in reference 1, we can find the reads that
### cover this position, then we can get the alignment position in reference 2 of
### these reads, and we can calculate the Pos2 in reference 2 that corresponds to
### the Pos1 in reference 1.


from .arg_parser import get_args_and_check_file
from .read_pos_resolver import get_reads_by_pos
from .ref_pos_resolver import get_new_pos_by_readname, get_new_pos_by_readname2, make_decision
from .input_output_handler import read_SNP_file, write_res1, open_res_file, close_res_file
import bri
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading


def process_snp(snp, args, bri_instance):
    ## get position of reads on reference1 
    res = get_reads_by_pos(args.bam1, args.ref1, snp)
    ## 注意判断是否为空列表,如果为空列表，就跳过
    if len(res) == 0:
        sys.stderr.write("No reads cover this position: %s\n" % snp)
        return

    ## get new position on reference2
    update_res = get_new_pos_by_readname(args.bam2, args.ref2, bri_instance, res)
    if len(update_res) == 0:
        sys.stderr.write("reads cant find in ref2 or all unmmaped: %s\n" % snp)
        return

    ## 做出判断决策
    final_res = make_decision(update_res)
    # with write_lock:
    #     write_res1(snp, final_res, outf)
    #     print(snp)
    #print(snp)
    return (snp, final_res)  # 修改此处，返回结果而不是写入文件

def snp_generator(snp_list):
    for snp in snp_list:
        yield snp

def main():
    ## 输出时间，用于计算运行时间
    start_time = time.time()
    sys.stderr.write("Start time: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time)))
    ## get args
    args = get_args_and_check_file()

    ## make the bri module ready, waiting for the read name that needed query
    bri_instance = bri.initialize_bri(args.bam2, args.bam2+".bri")

    ## 输出时间，用于计算运行时间
    sys.stderr.write("BRI module ready, time: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())))

    # outf = open_res_file(args.outprefix)
    # write_lock = threading.Lock()

    ## get snp position
    # snp_file = "/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/zz.CS1toCS2.1/06.small_batch_test/pipeline/test.txt"
    snp_list = read_SNP_file(args.snp)

    # 获取线程数
    num_threads = args.thread

    results = []  # 用于存储每个线程返回的结果

    # 初始化ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # 使用生成器
        tasks = {executor.submit(process_snp, snp, args, bri_instance): snp for snp in snp_generator(snp_list)}
        
        # 通过as_completed获取完成的任务
        for future in as_completed(tasks):
            try:
                result = future.result()
                if result:  # 如果结果不为空
                    results.append(result)
            except Exception as e:
                print(f"SNP {tasks[future]} generated an exception: {e}")


    # close_res_file(outf)
    with open_res_file(args.outprefix) as outf:
        for snp, final_res in results:
            write_res1(snp, final_res, outf)

    ## 输出时间，用于计算运行时间
    sys.stderr.write("End time: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())))
