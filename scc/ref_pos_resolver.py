import pandas as pd
import pysam
import bri
import tempfile, os

def get_new_reference_position(read, pos_in_read):
    cigar_tuples = read.cigartuples
    read_pos = 0
    ref_pos = read.reference_start

    # If the read is reverse complemented, adjust pos_in_read
    if read.is_reverse:
        pos_in_read = read.query_length - pos_in_read - 1

    for op, length in cigar_tuples:
        if op in [0, 7, 8]:  # M, =, X
            if read_pos <= pos_in_read < read_pos + length:
                return ref_pos + (pos_in_read - read_pos)   # 1-based
            read_pos += length
            ref_pos += length
        elif op == 1:  # I
            read_pos += length
        elif op == 2:  # D
            ref_pos += length
        elif op == 3:  # N
            ref_pos += length
        elif op == 4:  # soft clip
            read_pos += length
        elif op == 5:  # hard clip
            pass

    return None

## 根据read name 来获取新的位置
## @reture reads_dict: {(read_name, R1/R2): [read_pos, chr, new_ref_pos], ...}
def get_new_pos_by_readname(bam2_path, bri_instance, reads_name):
    ## 将其转为字典，read name 和 R1、R2 作为key，pos 作为value
    reads_dict = {(tuple(x[0:2])): [x[2],] for x in reads_name}

    ## 获取reads_name的第一个字段，为read name
    read_names_for_bri = list(set([x[0] for x in reads_name]))

    ## sam format text
    resulte = bri.query_by_readname(bam2_path, bri_instance, read_names_for_bri)
    # print(resulte)

    ## 将sam format text 转为 io.StringIO
    with tempfile.NamedTemporaryFile(delete=False) as temp:
        temp.write(resulte.encode())
        temp_path = temp.name

    # 使用 pysam 读取这个 "文件"
    with pysam.AlignmentFile(temp_path, "r") as samfile:
        # Iterate over the reads in the SAM file
        for read in samfile:
            read_name = read.query_name
            r1_or_r2 = 'R1' if read.is_read1 else 'R2'
            if (read_name, r1_or_r2) in reads_dict:
                pos = reads_dict[(read_name, r1_or_r2)][0]
                # Compute the chromosome and position in the new reference sequence
                chromosome = read.reference_name
                new_ref_pos = get_new_reference_position(read, pos)  # 1-based
                reads_dict[(read_name, r1_or_r2)].append(chromosome)
                reads_dict[(read_name, r1_or_r2)].append(new_ref_pos)

    os.remove(temp_path)


    # print(reads_dict)
    return reads_dict

## 根据多条reads 计算出来的新位置，来决定最终的位置
## @param reads_dicts: {(read_name, R1/R2): [read_pos, chr, new_ref_pos], ...}
## @return: [chr, new_ref_pos]
def make_decision(reads_dicts):
    ## 返回出现最多次数的chr 和 new_ref_pos的组合
    ## 如果出现次数相同，那么返回第一个
    ## 如果没有，那么返回None
    df = pd.DataFrame.from_dict(reads_dicts, orient='index', columns=['read_pos', 'chr', 'new_ref_pos'])
    df = df.groupby(['chr', 'new_ref_pos']).size().reset_index(name='counts')
    df = df.sort_values(by=['counts'], ascending=False)
    if df.shape[0] == 0:
        return None
    else:
        return [df.iloc[0]['chr'], df.iloc[0]['new_ref_pos']]


