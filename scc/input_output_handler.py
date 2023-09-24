
def read_SNP_file(snp_file):
    snp_list = []
    with open(snp_file,"r") as f:
        for line in f:
            snp_list.append(line.strip())

    return snp_list

def open_res_file(outprefix):
    f=  open(f"{outprefix}.txt","w")
    f.write("Ref1_chr\tRef1_pos\tRef1_base\tRef2_chr\tRef2_pos\tRef2_base\tcode\n")
    return f

def close_res_file(res_file):
    res_file.close()

def write_res1(snp, res, file_handle):
    ## snp: chr1B_part1:100
    ## res:['G', 'Chr1B_part1', 64796436, 'C', 1]
    snp_split = snp.split(":")
    file_handle.write(f"{snp_split[0]}\t{snp_split[1]}\t{res[0]}\t{res[1]}\t{res[2]}\t{res[3]}\t{res[4]}\n")
