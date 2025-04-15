import argparse
import sys,os

def get_args_and_check_file():
    cmdparser = argparse.ArgumentParser(
        description="SequenceCoordinateConverter (SCC) / fengcong@caas.cn" 
        )
    
    # 创建子命令分组
    subparsers = cmdparser.add_subparsers(dest='command', help='Commands')
    
    # 主程序转换坐标的参数
    convert_parser = subparsers.add_parser('convert', help='Convert coordinates')
    convert_parser.add_argument("-b","--bam1", dest="bam1",type=str, required=True,
        help="Bam/cram file of reference 1")
    convert_parser.add_argument("-B","--bam2", dest="bam2",type=str, required=True,
        help="Bam file of reference 2")
    convert_parser.add_argument("-r","--ref1", dest="ref1",type=str, required=True,
        help="Fasta file of reference 1")
    convert_parser.add_argument("-R","--ref2", dest="ref2",type=str, required=True,
        help="Fasta file of reference 2")
    convert_parser.add_argument("-s","--snp", dest="snp",type=str, required=True,
        help="SNP position file of reference 1, format: chr1A:100, one line one SNP")
    convert_parser.add_argument("-t","--thread", dest="thread",type=int, default=1,
        help="Thread number")
    convert_parser.add_argument("-o","--outprefix", dest="outprefix",type=str, required=True,
        help="Output prefix")
    convert_parser.add_argument("--use-server", dest="use_server", action="store_true",
        help="Try to use server mode if available")
    
    # 服务器模式的参数
    server_parser = subparsers.add_parser('server', help='Run in server mode')
    server_parser.add_argument("-B","--bam2", dest="bam2",type=str, required=True,
        help="Bam file of reference 2")
    server_parser.add_argument("-d","--daemonize", dest="daemonize", action="store_true",
        help="Run as daemon in background")
    
    # 停止服务器的参数
    stop_parser = subparsers.add_parser('stop-server', help='Stop server mode')
    
    args = cmdparser.parse_args()
    
    # 如果没有指定命令，默认使用转换坐标
    if args.command is None:
        args.command = 'convert'
        
    # 检查文件
    if args.command == 'convert':
        check_convert_files(args)
    elif args.command == 'server':
        check_server_files(args)
    
    return args

def check_convert_files(args):
    ## 检查每个文件是否存在，输出错误的话应该在 标准错误流
    ## 检查bam1是否是bam或者cram文件
    if not args.bam1.endswith(".cram") and not args.bam1.endswith(".bam") :
        sys.stderr.write("Error: reference1 must be bam or cram file!\n")
        sys.exit(1)

    ## 检查bam1是否存在
    if not os.path.exists(args.bam1):
        sys.stderr.write("Error: Bam file of reference 1 not exists!\n")
        sys.exit(1)
    
    ## 检查bam1的索引文件是否存在
    index_file1 = f"{args.bam1}.bai" if args.bam1.endswith(".bam") else f"{args.bam1}.crai"
    index_file2 = f"{args.bam1}.csi" if args.bam1.endswith(".bam") else f"{args.bam1}.crai"
    if not os.path.exists(index_file1) and not os.path.exists(index_file2):
        sys.stderr.write(f"Error: Index file {index_file1} does not exist!\n")
        sys.exit(1)

    ## 检查bam2 和其 索引 是否存在
    if not os.path.exists(args.bam2) or not os.path.exists(args.bam2+".bri"):
        sys.stderr.write("Error: Bam file or index file of reference 2 not exists!\n")
        sys.exit(1)
    
    ## 检查要转换的位置的文件是否存在
    if not os.path.exists(args.snp):
        sys.stderr.write("Error: SNP position file not exists!\n")
        sys.exit(1)

    ## 如果提供了ref1 和 ref2 就检查其文件是否存在, 还有其索引文件
    if not os.path.exists(args.ref1) or not os.path.exists(args.ref1+".fai"):
        sys.stderr.write("Error: Fasta file and index of reference 1 not exists!\n")
        sys.exit(1)
    if not os.path.exists(args.ref2) or not os.path.exists(args.ref2+".fai"):
        sys.stderr.write("Error: Fasta file and index of reference 2 not exists!\n")
        sys.exit(1)

def check_server_files(args):
    ## 检查bam2 和其 索引 是否存在
    if not os.path.exists(args.bam2) or not os.path.exists(args.bam2+".bri"):
        sys.stderr.write("Error: Bam file or index file of reference 2 not exists!\n")
        sys.exit(1)

