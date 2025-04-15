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
import os
import signal
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import atexit

# 全局变量，用于服务器模式
shm_id = None
shm_addr = None
bri_instance_ptr = None
client_ctx_ptr = None  # 添加客户端上下文指针

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
    return (snp, final_res)  # 修改此处，返回结果而不是写入文件

def snp_generator(snp_list):
    for snp in snp_list:
        yield snp

def daemonize():
    """守护进程化"""
    try:
        pid = os.fork()
        if pid > 0:
            # 退出父进程
            sys.exit(0)
    except OSError as e:
        sys.stderr.write(f"Fork #1 failed: {e}\n")
        sys.exit(1)
    
    # 分离父目录
    os.chdir("/")
    os.setsid()
    os.umask(0)
    
    # 第二次fork
    try:
        pid = os.fork()
        if pid > 0:
            # 退出第一个子进程
            sys.exit(0)
    except OSError as e:
        sys.stderr.write(f"Fork #2 failed: {e}\n")
        sys.exit(1)
    
    # 重定向标准文件描述符
    sys.stdout.flush()
    sys.stderr.flush()
    
    si = open(os.devnull, 'r')
    so = open(os.devnull, 'a+')
    se = open(os.devnull, 'a+')
    
    os.dup2(si.fileno(), sys.stdin.fileno())
    os.dup2(so.fileno(), sys.stdout.fileno())
    os.dup2(se.fileno(), sys.stderr.fileno())

def cleanup_server():
    """清理服务器资源"""
    try:
        # 如果是因为异常或者信号导致退出，可能会执行多次cleanup
        # 忽略错误信息避免在日志中出现多余错误
        bri.stop_server()
        sys.stderr.write("Server stopped, shared memory cleaned up.\n")
    except Exception as e:
        # 只在调试模式下或在非正常退出时打印详细错误
        if str(e) != "Server not running or already stopped":
            sys.stderr.write(f"Error stopping server: {e}\n")

def disconnect_from_shared_memory():
    """从共享内存断开连接"""
    global shm_addr, client_ctx_ptr
    if shm_addr is not None:
        try:
            sys.stderr.write(f"Disconnecting from shared memory at address {shm_addr}\n")
            # 传递客户端上下文以便释放
            bri.disconnect_from_server(shm_addr, client_ctx_ptr)
            shm_addr = None
            client_ctx_ptr = None
        except Exception as e:
            sys.stderr.write(f"Error disconnecting from shared memory: {e}\n")
            # 即使出错，也将指针设为None以避免重复尝试断开连接
            shm_addr = None
            client_ctx_ptr = None

def signal_handler(sig, frame):
    """处理信号，确保干净退出"""
    sys.stderr.write("Received signal to terminate.\n")
    cleanup_server()
    sys.exit(0)

def run_server(args):
    """运行服务器模式"""
    global shm_id
    
    sys.stderr.write("Starting server mode...\n")
    
    # 先检查是否已有服务器在运行
    key_path = "/tmp/scc_shm_key"
    size_path = "/tmp/scc_shm_size"
    if os.path.exists(key_path) or os.path.exists(size_path):
        sys.stderr.write("Warning: Previous server session detected. Cleaning up...\n")
        try:
            cleanup_server()
        except Exception:
            # 如果清理失败，可能需要手动删除文件
            if os.path.exists(key_path):
                os.unlink(key_path)
            if os.path.exists(size_path):
                os.unlink(size_path)
    
    # 如果需要，变成守护进程
    if args.daemonize:
        sys.stderr.write("Daemonizing...\n")
        daemonize()
    
    # 注册信号处理程序
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    # 注册退出清理
    atexit.register(cleanup_server)
    
    try:
        # 初始化BRI并创建共享内存
        sys.stderr.write(f"Loading index from {args.bam2}...\n")
        shm_id = bri.create_server(args.bam2, args.bam2+".bri")
        
        sys.stderr.write(f"Server started with shared memory ID: {shm_id}\n")
        sys.stderr.write("Press Ctrl+C to stop server.\n")
        
        # 保持程序运行
        while True:
            time.sleep(60)
    
    except Exception as e:
        sys.stderr.write(f"Error in server mode: {e}\n")
        cleanup_server()
        sys.exit(1)

def stop_server():
    """停止服务器模式"""
    try:
        # 尝试获取key文件路径
        key_path = "/tmp/scc_shm_key"
        if not os.path.exists(key_path):
            sys.stderr.write("Server is not running.\n")
            return
        
        # 发送终止信号给服务器进程
        # 首先尝试找到服务器进程
        server_pid = None
        for pid in [int(pid) for pid in os.listdir('/proc') if pid.isdigit()]:
            try:
                cmdline = open(f'/proc/{pid}/cmdline', 'rb').read().decode().replace('\0', ' ')
                if 'SCC server' in cmdline:
                    server_pid = pid
                    break
            except (IOError, ProcessLookupError):
                continue
        
        if server_pid:
            os.kill(server_pid, signal.SIGTERM)
            sys.stderr.write(f"Sent termination signal to server process (PID: {server_pid}).\n")
            # 给进程一些时间来执行cleanup
            time.sleep(1)
        else:
            sys.stderr.write("No running server process found. Cleaning up manually.\n")
            cleanup_server()
    
    except Exception as e:
        # 如果是因为服务器已经停止而导致的错误，不要显示错误
        if "Server not running" not in str(e):
            sys.stderr.write(f"Error stopping server: {e}\n")
        sys.exit(1)

def run_convert(args):
    """运行坐标转换"""
    global shm_addr, bri_instance_ptr, client_ctx_ptr
    
    start_time = time.time()
    sys.stderr.write("Start time: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time)))
    
    bri_instance = None
    
    # 如果启用了服务器模式，尝试连接
    if hasattr(args, 'use_server') and args.use_server:
        sys.stderr.write("Attempting to connect to server...\n")
        try:
            result = bri.connect_to_server(args.bam2)
            if result is not None and len(result) >= 3:
                shm_addr, bri_instance_ptr, client_ctx_ptr = result
                bri_instance = bri_instance_ptr
                # 注册退出时断开连接
                atexit.register(disconnect_from_shared_memory)
                sys.stderr.write(f"Connected to server successfully. SHM address: {shm_addr}, BRI pointer: {bri_instance_ptr}, Context pointer: {client_ctx_ptr}\n")
                
                # 安全检查 - 验证BRI结构
                try:
                    # 检查BRI指针是否还是有效的
                    sys.stderr.write(f"Validating BRI structure at {bri_instance}...\n")
                    
                    # 增强型测试 - 检查是否能够执行一次简单查询
                    test_read_name = "test_dummy_read"
                    sys.stderr.write(f"Testing with read name: {test_read_name}\n")
                    test_result = bri.query_by_readname(args.bam2, bri_instance, [test_read_name])
                    sys.stderr.write(f"BRI validation test completed successfully\n")
                except Exception as e:
                    sys.stderr.write(f"ERROR: BRI validation test failed: {str(e)}\n")
                    sys.stderr.write("Falling back to direct index loading\n")
                    bri_instance = None
                    disconnect_from_shared_memory()
            else:
                sys.stderr.write("Failed to connect to server. Will load index directly.\n")
        except Exception as e:
            sys.stderr.write(f"Error connecting to server: {str(e)}\n")
            sys.stderr.write("Will load index directly.\n")
            # 确保断开共享内存连接，避免内存泄露
            if shm_addr is not None:
                try:
                    disconnect_from_shared_memory()
                except:
                    pass
    
    # 如果没有连接到服务器，或者不使用服务器模式，则常规加载索引
    if bri_instance is None:
        sys.stderr.write("Loading index directly...\n")
        try:
            bri_instance = bri.initialize_bri(args.bam2, args.bam2+".bri")
            sys.stderr.write(f"Index loaded directly, BRI pointer: {bri_instance}\n")
        except Exception as e:
            sys.stderr.write(f"Error loading index: {str(e)}\n")
            sys.exit(1)
    
    sys.stderr.write("BRI module ready, time: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())))
    
    # 获取SNP位置
    snp_list = read_SNP_file(args.snp)
    sys.stderr.write(f"Loaded {len(snp_list)} SNP positions\n")
    
    # 获取线程数
    num_threads = args.thread
    sys.stderr.write(f"Using {num_threads} threads\n")
    
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
                sys.stderr.write(f"SNP {tasks[future]} generated an exception: {str(e)}\n")
                import traceback
                traceback.print_exc(file=sys.stderr)
    
    sys.stderr.write(f"Successfully processed {len(results)} SNP positions\n")
    
    # 写入结果
    outf = open_res_file(args.outprefix)
    for snp, final_res in results:
        write_res1(snp, final_res, outf)
    close_res_file(outf)
    sys.stderr.write(f"Results written to {args.outprefix}.txt\n")
    
    # 如果使用了共享内存，断开连接
    if shm_addr is not None:
        sys.stderr.write("Disconnecting from shared memory\n")
        disconnect_from_shared_memory()
    
    # 输出时间，用于计算运行时间
    end_time = time.time()
    sys.stderr.write("End time: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time)))
    sys.stderr.write(f"Total execution time: {end_time - start_time:.2f} seconds\n")

def main():
    # 获取参数
    args = get_args_and_check_file()
    
    # 根据命令运行相应功能
    if args.command == 'convert':
        run_convert(args)
    elif args.command == 'server':
        run_server(args)
    elif args.command == 'stop-server':
        stop_server()
    else:
        sys.stderr.write(f"Unknown command: {args.command}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
