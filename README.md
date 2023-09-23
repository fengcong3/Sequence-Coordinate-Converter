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

b. Sort the alignment file for `reference1` (`alignment1`) based on genomic coordinates and create an index using `samtools`.

```bash
samtools sort -o alignment1_sorted.bam -O bam alignment1.sam
samtools index alignment1_sorted.bam
```

c. Sort the alignment file for `reference2` (`alignment2`) based on the read name.

```bash
samtools sort -n -o alignment2_namesort.bam -O bam alignment2.sam
```

d. Use `bri` to index the sorted `alignment2`.

```bash
bri index alignment2_namesort.bam
```

e. Use `SCC` for position conversion.

```bash
python main.py --input alignment1_sorted.bam --cross alignment2_namesort.bam --output converted_coordinates.txt
```

### 2. Coordinate Conversion Between Different Species or Closely Related Species

a. Choose the reads of either `reference1` or `reference2` and align them to both references using `bwa mem`.

```bash
bwa mem reference1.fa chosen_reads.fq > alignment1.sam
bwa mem reference2.fa chosen_reads.fq > alignment2.sam
```

... (Rest of the steps would be similar to the first scenario)

---

## Discussion

For the second scenario, it's not entirely clear whether one should use reads from `reference1` or `reference2` for the conversion, or perhaps both. This is an area that requires further investigation. The choice might depend on specific use-cases or the quality and coverage of the available reads.

By selecting reads from `reference1`, you ensure that the coordinates are more directly comparable to `reference2`. On the other hand, using reads from `reference2` might provide more context for the alignment, thereby potentially improving the accuracy of the conversion. Utilizing reads from both references might offer a comprehensive approach but could complicate the analysis.

A side-by-side comparison or hybrid approach employing both sets of reads could offer more insights into this issue.