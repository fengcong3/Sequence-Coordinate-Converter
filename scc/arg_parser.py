import argparse
import sys,os

def get_args_and_check_file():
    cmdparser = argparse.ArgumentParser(
        description="SequenceCoordinateConverter (SCC) / fengcong@caas.cn" 
        )
    cmdparser.add_argument("-b","--bam1", dest="bam1",type=str, required=True,
        help="Bam file of reference 1")
    cmdparser.add_argument("-B","--bam2", dest="bam2",type=str, required=True,
        help="Bam/cram file of reference 2")
    cmdparser.add_argument("-r","--ref1", dest="ref1",type=str, required=False,
        help="Fasta file of reference 1")
    cmdparser.add_argument("-R","--ref2", dest="ref2",type=str, required=False,
        help="Fasta file of reference 2")
    cmdparser.add_argument("-s","--snp", dest="snp",type=str, required=True,
        help="SNP position file of reference 1, format: chr1A:100, one line one SNP")
    ## thread
    cmdparser.add_argument("-t","--thread", dest="thread",type=int, default=1,
        help="Thread number")
    cmdparser.add_argument("-o","--outprefix", dest="outprfix",type=str, required=True,
        help="Output prefix")
    
    args = cmdparser.parse_args()

    ## 检查每个文件是否存在，输出错误的话应该在 标准错误流
    ## 检查bam1是否是bam或者cram文件
    if not args.bam1.endswith(".cram") and args.bam1.endswith(".bam") :
        sys.stderr.write("Error: reference1 must be bam or cram file!\n")
        sys.exit(1)

    ## 检查bam1是否存在
    if not os.path.exists(args.bam1):
        sys.stderr.write("Error: Bam file of reference 1 not exists!\n")
        sys.exit(1)
    
    ## 检查bam1的索引文件是否存在
    index_file = f"{args.bam1}.bai" if args.bam1.endswith(".bam") else f"{args.bam1}.crai"
    if not os.path.exists(index_file):
        sys.stderr.write(f"Error: Index file {index_file} does not exist!\n")
        sys.exit(1)

    ## 检查bam2 和其 索引 是否存在
    if not os.path.exists(args.bam2) or not os.path.exists(args.bam2+".bri"):
        sys.stderr.write("Error: Bam file or index file of reference 2 not exists!\n")
        sys.exit(1)
    
    ## 检查要转换的位置的文件是否存在
    if not os.path.exists(args.snp):
        sys.stderr.write("Error: SNP position file not exists!\n")
        sys.exit(1)

    ## 如果bam1 是cram文件，那么必须提供ref1
    if args.bam1.endswith(".cram") and not args.ref1:
        sys.stderr.write("Error: If bam file is cram, you must provide reference 1!\n")
        sys.exit(1)
    
    ## 如果提供了ref1 和 ref2 就检查其文件是否存在
    if args.ref1 and not os.path.exists(args.ref1):
        sys.stderr.write("Error: Fasta file of reference 1 not exists!\n")
        sys.exit(1)
    if args.ref2 and not os.path.exists(args.ref2):
        sys.stderr.write("Error: Fasta file of reference 2 not exists!\n")
        sys.exit(1)

    return args
    