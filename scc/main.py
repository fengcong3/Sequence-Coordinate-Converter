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
from .ref_pos_resolver import get_new_pos_by_readname, make_decision
from .input_output_handler import read_SNP_file, write_res1, open_res_file, close_res_file
import bri
import sys

def main():
    ## get args
    args = get_args_and_check_file()

    ## make the bri module ready, waiting for the read name that needed query
    # bam2="/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/zz.CS1toCS2.1//01.cram_of_CS2/CS/CS2.merge.sortn.bam"
    # ref2="/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/zz.CS1toCS2.1/01.cram_of_CS2/CS/CS_V2.1.cut.norm.fa"
    bri_instance = bri.initialize_bri(args.bam2, args.bam2+".bri")

    outf = open_res_file(args.outprefix)

    ## get snp position
    # snp_file = "/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/zz.CS1toCS2.1/06.small_batch_test/pipeline/test.txt"
    snp_list = read_SNP_file(args.snp)

    ## traverse each SNP
    # bam_file_path = "/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/zz.CS1toCS2.1/00.cram_of_CS1/merge.rmdup.Chinese_Spring.cram"
    # ref1="/public/agis/chengshifeng_group/fengcong/WGRS/graduation_project/00.var_genome/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta"
    for snp in snp_list:
        ## get position of reads on reference1 
        res=get_reads_by_pos(args.bam1, args.ref1, snp)
        ## 注意判断是否为空列表,如果为空列表，就跳过
        if len(res) == 0:
            sys.stderr.write("No reads cover this position: %s\n" % snp)
            continue

        ## get new position on reference2
        update_res = get_new_pos_by_readname(args.bam2, args.ref2, bri_instance, res)
        if len(update_res) == 0:
            sys.stderr.write("reads cant find in ref2 or all unmmaped: %s\n" % snp)
            continue

        ## 做出判断决策
        final_res=make_decision(update_res)
        write_res1(snp, final_res, outf)

    close_res_file(outf)