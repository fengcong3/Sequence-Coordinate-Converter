import pysam
# import sys


def get_read_position_in_reference1(read, ref_pos):
    # start position of the read in reference
    start = read.pos
    # end position of the read in reference
    end = read.pos + read.alen - 1

    if ref_pos < start or ref_pos > end:
        return None

    ref_cursor = start
    read_cursor = 0 if not read.is_reverse else read.query_length - 1

    for op, length in read.cigar:
        # 'M' operation
        if op == 0:
            if ref_cursor <= ref_pos <= ref_cursor + length - 1:
                offset = ref_pos - ref_cursor
                return read_cursor + offset if not read.is_reverse else read_cursor - offset
            ref_cursor += length
            read_cursor += length if not read.is_reverse else -length

        # 'I' operation
        elif op == 1:
            read_cursor += length if not read.is_reverse else -length

        # 'D' or 'N' operation
        elif op in [2, 3]:
            if ref_cursor <= ref_pos <= ref_cursor + length - 1:
                return None
            ref_cursor += length

        # 'S' or 'H' operation
        elif op in [4, 5]:
            read_cursor += length if not read.is_reverse else -length

    return None

## 提取指定位置的reads，如果是cram 可能需要reference
## 然后根据 这些reads 和 参考基因组的位置，来计算在read 上的位置
## @param bam_file_path: bam file path
## @param ref1: reference 1
## @param pos:  chr1A_part1:100
## @return: [[read_name, R1/R2, read_pos], ...]
def get_reads_by_pos(bam_file_path, ref1, pos):
    # 根据文件类型选择打开模式
    mode = "rb" if bam_file_path.endswith( ".bam" )  else "rc"

    # 可选：对于CRAM，指定参考基因组
    reference_file = ref1 if bam_file_path.endswith( ".cram" )  else None

    # 处理pos
    chr, pos = pos.split(":")
    pos = int(pos)

    # 定义返回值
    reads = []

    # 打开文件
    with pysam.AlignmentFile(bam_file_path, mode, reference_filename=reference_file) as infile:
        # 遍历覆盖特定位置（例如chr1A:100-100）的reads
        # 注意：pysam的坐标是0-based，而bam文件里面是1-based
        for read in infile.fetch(chr, pos-1, pos):
            if not read.is_unmapped:
                read_pos = get_read_position_in_reference1(read, pos)
                if read_pos is not None:
                    r1_or_r2 = 'R1' if read.is_read1 else 'R2'
                    # print(f'{read.query_name},{r1_or_r2},{read_pos}')
                    reads.append([read.query_name, r1_or_r2, read_pos])

    return reads
            
