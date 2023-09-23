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

import pandas as pd

from .arg_parser import get_args_and_check_file
from .read_pos_resolver import get_reads_by_pos
import bri
import io

def main():
    ## get args
    # args = get_args_and_check_file()

    ## read_pos_resolver.py
    bam_file_path = "/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/zz.CS1toCS2.1/00.cram_of_CS1/merge.rmdup.Chinese_Spring.cram"
    ref1="/public/agis/chengshifeng_group/fengcong/WGRS/graduation_project/00.var_genome/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta"
    pos = "chr1B_part1:882559"
    res=get_reads_by_pos(bam_file_path, ref1, pos)
    print(res)

    ## make the bri module ready, waiting for the read name that needed query
    # bri_instance = bri.initialize_bri(args.bam2, args.bam2+".bri")



    # 假设read_names是一个Python列表
    # input_bam_path= bam_file_path
    # read_names = ["SRR5893652.16936768", "SRR5893652.185558604" ]
    # output_sam_path = "/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/zz.CS1toCS2.1/06.small_batch_test/pipeline/zz.scctest.sam"
    # resulte = bri.query_by_readname(input_bam_path, bri_instance, read_names)
    # print(resulte)

    # sam_file = io.StringIO(resulte)

    # # 使用 pysam 读取这个 "文件"
    # with pysam.AlignmentFile(sam_file, "r") as samfile:
    #     for read in samfile:
    #         print(read)
