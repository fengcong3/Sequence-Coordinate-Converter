# SequenceCoordinateConverter (SCC)
A tool for converting genomic coordinates between different versions of the same species or between closely related species.

## Method
```
--------|---------reference 1----------------- chr1A
    ----|-read1----->
    <---|--read2-----
      --|---read3----->

## "|" represents the pos need to be converted

-----------------reference 2------/----------- chr1A
                              ----/-read1----->
                              <---/--read2-----
                                --/---read3----->

## "/" represents the pos has been converted
```
When we have a position (Pos1) in reference 1, we can find the reads that cover this position, then we can get the alignment position in reference 2 of these reads, and we can calculate the Pos2 in reference 2 that corresponds to the Pos1 in reference 1.


## Usage

### 1. Coordinate Conversion for the Same Species with Different Assembly Versions

a. Align the second-generation reads of the species to both `reference1` and `reference2` using `bwa mem`.

```bash
bwa mem reference1.fa reads.fq > alignment1.sam
bwa mem reference2.fa reads.fq > alignment2.sam
```

b. Sort the alignment file for `reference1` based on genomic coordinates and create an index using `samtools`.

```bash
samtools sort -o alignment1_sorted.bam -O bam alignment1.sam
samtools index alignment1_sorted.bam
```

c. Sort the alignment file for `reference2` based on the read name.

```bash
samtools sort -n -o alignment2_namesort.bam -O bam alignment2.sam
```

d. Use `bri` to index the sorted `alignment2`.

```bash
bri index alignment2_namesort.bam
```

e. Use `SCC` for position conversion.

```bash
$ SCC -h
usage: SCC [-h] -b BAM1 -B BAM2 -r REF1 -R REF2 -s SNP -o OUTPREFIX

SequenceCoordinateConverter (SCC) / fengcong@caas.cn

optional arguments:
  -h, --help            show this help message and exit
  -b BAM1, --bam1 BAM1  Bam file of reference 1
  -B BAM2, --bam2 BAM2  Bam/cram file of reference 2
  -r REF1, --ref1 REF1  Fasta file of reference 1
  -R REF2, --ref2 REF2  Fasta file of reference 2
  -s SNP, --snp SNP     SNP position file of reference 1, format: chr1A:100,
                        one line one SNP
  -o OUTPREFIX, --outprefix OUTPREFIX
                        Output prefix
```

### 2. Coordinate Conversion Between Different Species or Closely Related Species

a. Choose the reads of either `reference1` or `reference2` and align them to both references using `bwa mem`.

```bash
bwa mem reference1.fa chosen_reads.fq > alignment1.sam
bwa mem reference2.fa chosen_reads.fq > alignment2.sam
```

... (Rest of the steps would be similar to the first scenario)

## Results

Example result:
```
Ref1_chr	Ref1_pos	Ref1_base	Ref2_chr	Ref2_pos	Ref2_base	code
chr1B_part1	882559	A	Chr1B_part1	875080	A	0
chr1B_part1	58917715	G	Chr1B_part1	64796436	C	1
```

### Explanation of the `code` Field

The `code` field is used to describe the relationship between `Ref1` and `Ref2` at the corresponding positions:

- `0`: No change in assembly between `Ref1` and `Ref2` at this position.
- `1`: Change in assembly orientation between `Ref1` and `Ref2` at this position.
- `2`: Difference in assembly between `Ref1` and `Ref2` at this position, but the orientation remains unchanged.
- `3`: Difference in assembly between `Ref1` and `Ref2` at this position, and the orientation has also changed.




## Discussion

For the second scenario, it's not entirely clear whether one should use reads from `reference1` or `reference2` for the conversion, or perhaps both. This is an area that requires further investigation.

A side-by-side comparison or hybrid approach employing both sets of reads could offer more insights into this issue.