import pandas as pd
import pysam
import bri
import tempfile, os
import re

# 解析CIGAR字符串为一个元组列表
def parse_cigar(cigar_str):
    ops = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
    op_tuples = []
    lengths = re.split("[MIDNSHP=X]", cigar_str)
    for i, op in enumerate(re.findall("[MIDNSHP=X]", cigar_str)):
        op_tuples.append((ops.index(op), int(lengths[i])))
    return op_tuples

# 根据CIGAR元组列表，计算新的参考位置
def get_new_reference_position2(cigar_tuples, pos_in_read, reference_start, is_reverse, query_length):
    read_pos = 0
    ref_pos = reference_start

    # 如果是反向互补，调整 pos_in_read
    if is_reverse:
        pos_in_read = query_length - pos_in_read + 1

    for op, length in cigar_tuples:
        if op in [0, 7, 8]:  # M, =, X
            if read_pos <= pos_in_read < read_pos + length:
                return ref_pos + (pos_in_read - read_pos)  # 1-based
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
    return None

## 根据read name 来获取新的位置
## @reture reads_dict: {(read_name, R1/R2): [read_pos,ref1_base,read_base, chr, new_ref_pos,ref2_base, read_base], ...}
def get_new_pos_by_readname2(bam2_path, ref2, bri_instance, reads_name):
    ## 将其转为字典，read name 和 R1、R2 作为key，pos 作为value
    reads_dict = {(tuple(x[0:2])): x[2:] for x in reads_name}

    ## 获取reads_name的第一个字段，为read name
    read_names_for_bri = list(set([x[0] for x in reads_name]))

    ## 打开参考基因组文件
    fasta = pysam.FastaFile(ref2)
    

    ## sam format text
    resulte = bri.query_by_readname(bam2_path, bri_instance, read_names_for_bri)
    # print(resulte)

    # 使用 Python 的 str.split 将结果分割为行
    sam_lines = resulte.strip().split('\n')

    for line in sam_lines:
        if line.startswith('@'):
            continue
        
        fields = line.split('\t')

        # 提取所需字段
        read_name = fields[0]
        flag = int(fields[1])
        reference_name = fields[2]
        reference_start = int(fields[3])  # 1-based
        mapping_quality = int(fields[4])
        cigar_str = fields[5]
        sequence = fields[9]
        query_length = len(sequence)

        is_reverse = bool(flag & 16)
        is_read1 = bool(flag & 64)
        is_supplementary = bool(flag & 2048)

        r1_or_r2 = 'R1' if is_read1 else 'R2'

        ## 如果比对质量太低，或者是supplementary alignment，就跳过
        if mapping_quality < 1 or is_supplementary or "SA" in line:
            continue

        if (read_name, r1_or_r2) in reads_dict:
            pos = reads_dict[(read_name, r1_or_r2)][0]
            # Compute the chromosome and position in the new reference sequence
            chromosome = reference_name

            # 解析CIGAR
            cigar_tuples = parse_cigar(cigar_str)

            # 计算新位置
            new_ref_pos = get_new_reference_position2(cigar_tuples, pos, reference_start, is_reverse, len(sequence))

            ## 会出现supplementary alignment 的情况，这种情况下，
            ## 会有两条record，一条是主要的，一条是supplementary的
            ## 那么有一条 就不会有位置？
            if new_ref_pos is None: ## 如果没找到位置，就跳过
                continue
            # print(read_name, r1_or_r2)
            # print(reads_dict[(read_name, r1_or_r2)])
            # print(pos, new_ref_pos)
            # print(read.query_sequence, read.query_length)
            # print(pos-1, read.query_length - pos)
            # 获取某个染色体(chr1)上第1000个位置（0-based）的碱基
            ref_base_at_position = fasta.fetch(chromosome, new_ref_pos-1, new_ref_pos).upper()
            read_base_at_position = sequence[pos-1].upper() if not is_reverse else sequence[query_length - pos ].upper()
            reads_dict[(read_name, r1_or_r2)].append(chromosome)
            reads_dict[(read_name, r1_or_r2)].append(new_ref_pos)
            reads_dict[(read_name, r1_or_r2)].append(ref_base_at_position)
            reads_dict[(read_name, r1_or_r2)].append(read_base_at_position)


    # 删除没有更新的键
    keys_to_remove = [key for key, value in reads_dict.items() if len(value) < 7]  # 假设一个完整的value应该包含8个元素
    for key in keys_to_remove:
        del reads_dict[key]
    # print(reads_dict)
    return reads_dict

def get_new_reference_position(read, pos_in_read):
    cigar_tuples = read.cigartuples
    read_pos = 0
    ref_pos = read.reference_start

    # If the read is reverse complemented, adjust pos_in_read
    if read.is_reverse:
        pos_in_read = read.query_length - pos_in_read +1

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
## @reture reads_dict: {(read_name, R1/R2): [read_pos,ref1_base,read_base, chr, new_ref_pos,ref2_base, read_base], ...}
def get_new_pos_by_readname(bam2_path, ref2, bri_instance, reads_name):
    ## 将其转为字典，read name 和 R1、R2 作为key，pos 作为value
    reads_dict = {(tuple(x[0:2])): x[2:] for x in reads_name}

    ## 获取reads_name的第一个字段，为read name
    read_names_for_bri = list(set([x[0] for x in reads_name]))

    ## 打开参考基因组文件
    fasta = pysam.FastaFile(ref2)
    

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

            ## 如果比对质量太低，或者是supplementary alignment，就跳过
            if read.mapping_quality < 1 or read.is_supplementary or "SA" in dict(read.tags):
                continue

            if (read_name, r1_or_r2) in reads_dict:
                pos = reads_dict[(read_name, r1_or_r2)][0]
                # Compute the chromosome and position in the new reference sequence
                chromosome = read.reference_name
                new_ref_pos = get_new_reference_position(read, pos)  # pos:1-based

                ## 会出现supplementary alignment 的情况，这种情况下，
                ## 会有两条record，一条是主要的，一条是supplementary的
                ## 那么有一条 就不会有位置？
                if new_ref_pos is None: ## 如果没找到位置，就跳过
                    continue
                # print(read_name, r1_or_r2)
                # print(reads_dict[(read_name, r1_or_r2)])
                # print(pos, new_ref_pos)
                # print(read.query_sequence, read.query_length)
                # print(pos-1, read.query_length - pos)
                # 获取某个染色体(chr1)上第1000个位置（0-based）的碱基
                ref_base_at_position = fasta.fetch(chromosome, new_ref_pos-1, new_ref_pos).upper()
                read_base_at_position = read.query_sequence[pos-1].upper() if not read.is_reverse else read.query_sequence[read.query_length - pos ].upper()
                reads_dict[(read_name, r1_or_r2)].append(chromosome)
                reads_dict[(read_name, r1_or_r2)].append(new_ref_pos)
                reads_dict[(read_name, r1_or_r2)].append(ref_base_at_position)
                reads_dict[(read_name, r1_or_r2)].append(read_base_at_position)

    os.remove(temp_path)
    del(resulte)
    # 删除没有更新的键
    keys_to_remove = [key for key, value in reads_dict.items() if len(value) < 7]  # 假设一个完整的value应该包含8个元素
    for key in keys_to_remove:
        del reads_dict[key]
    # print(reads_dict)
    return reads_dict

## 根据多条reads 计算出来的新位置，来决定最终的位置
## @param reads_dicts: {(read_name, R1/R2): [read_pos, ref1_base, read_base1, chr, new_ref_pos,ref2_base, read_base2], ...}
## @return:  [ref1_base，chr，new_ref_pos，ref2_base，标识符 ]  
## 其中标识符 0代表碱基没有发生改变，1代表两个参考基因组组装方向发生改变，变为互补碱基，2代表两个基因组组装结果不一致;3代表反方向的不一致
## warning: 仅用了最热点的一个位置的一条read 做判断，不清楚有啥后果
def make_decision(reads_dicts):
    df = pd.DataFrame.from_dict(reads_dicts, orient='index', 
                                columns=['read_pos', 'ref1_base', 'read_base1', 
                                         'chr', 'new_ref_pos', 'ref2_base', 'read_base2'])
                                         
    # 对['chr', 'new_ref_pos']进行分组，并计算每组的大小
    grouped = df.groupby(['chr', 'new_ref_pos'])
    
    # 获取出现次数最多的['chr', 'new_ref_pos']组合
    most_common_group = grouped.size().idxmax()
    filtered_df = grouped.get_group(most_common_group)
    
    # 初始化结果列表
    result = [filtered_df['ref1_base'].iloc[0], most_common_group[0], most_common_group[1], filtered_df['ref2_base'].iloc[0]]
    
    # 判断碱基状态的改变
    ref1_base = filtered_df['ref1_base'].iloc[0]
    read_base1 = filtered_df['read_base1'].iloc[0]
    ref2_base = filtered_df['ref2_base'].iloc[0]
    read_base2 = filtered_df['read_base2'].iloc[0]

    # 标识符初始化
    identifier = 0

    # 判断是否变为互补碱基
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if complement_dict.get(ref1_base) == ref2_base and complement_dict.get(read_base1) == read_base2:
        identifier = 1
    # 基因组反向且不一致
    elif ref1_base != ref2_base and complement_dict.get(read_base1) == read_base2:
        identifier = 3
    # 判断是否两个基因组组装结果不一致
    elif ref1_base != ref2_base or read_base1 != read_base2:
        identifier = 2

    result.append(identifier)
    
    return result





