# Title: NGS data coverage analysis for NA12878
# Author: Ying D.
# Date: 10-24-2024

## Project description 

Get an NGS BAM file generated from whole genome sequencing (WGS) data aligned against the Genome in a Bottle (GIAB, NA12878/HG001) reference sample. The BAM file contains aligned reads from the sequencing process. Each read represents a fragment of DNA from the genome. The coverage at a specific genomic position is defined as the number of reads that overlap that position.
Please calculate the coverage(sequencing depth) for the given NA12878 BAM (https://www.internationalgenome.org/data-portal/sample/NA12878) file aligned against GRCH38 reference sequence. 

## Analytical plan design

This analysis aims to assess read coverage across the entire genome (HG38) for the sample NA12878. My proposed analysis plan has two main aims.

The first aim is to examine the read coverage at various genomic positions, and determine the average coverage across the genome with the non-zero coverage. To accomplish this, we need to generate a coverage file for various genomic positions and calculate the following: the number of unique genomic positions with read coverage, the total number of reads that cover these positions, and the total length of the combined genomic positions. We can then use these values to calculate the average coverage.

The second aim is to analyze the distribution of read coverage across the genome, which will help to assess whether read coverage is evenly distributed throughout the genome or if there are any biases present. To accomplish this, we need to generate a histogram of coverage file to determine the fraction of bases on the chromosomes or the entire genome, along with their corresponding depth of coverage

The detailed plan aims and objectives are listed below.

### Aims: 
#### Aim 1: Examine the read coverage at various genomic location and determine the average read coverage depth across the genome with the non-zero coverage.
***Objectives:*** 
* Generate a coverage file that displays the non-zero coverage for each position in the genome (base pairs with 0 coverage are excluded).
* Count the total number of genomic positions with non-zero coverage.
* Calculate the total number of reads by summing all the reads across genomic positions with the non-zero coverage.
* Calculate the average read coverage per read-covered genomic position.
   ```
   Formula: Coverage per sequenced position = Total_read_number / Total_postion_number
   ```
* Calculate the total length of combined genomic postions with non-zero coverage (sequenced_genome_length) in base pairs (bp) by summing the length of all genomic positions with the non-zero coverage.
* Calculate the average read coverage per base pair on the sequenced genome.
   ```
   Formula: Coverage per base pair = Total_read_number / sequenced_genome_length
   ```
* Calculate the proportion of the whole genome that has read coverage, given that the length of the HG38 (GRCh38) genome is approximately 3,200,000,000 base pairs (3.2 billion base pairs).
    ```
    Formula: sequenced_genome_proportion = sequenced_genome_length / 3.2*10^9​
    ```
#### Aim 2: Analyze the distribution of read coverage across the genome
***Objectives:***
* Generate a histogram of coverage data that shows the depth of coverage across chromosomes or the entire genome, along with the fraction of bases on chromosomes (or the entire genome) with corresponding depths.
* Determine the distribution of read coverage on chromosomes or the entire genome.

## Dataset: 
NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram

source location: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram
> NOTE: 
> There are two suggested datasets for this project: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bam.bas and ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram. 
The first file is not a BAM file; however, there is a CRAM file with the same name in the same directory.  The second dataset is too large (15 GB), making it difficult to handle on my local machine. Therefore, this project used the .cram file from the first directory for analysis.

## Approach & Algorithms: 
#### 1. Download the data

```
curl 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram' > 'Data/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram'

curl 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram.crai' > 'Data/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram.crai'
```

#### 2. Generate a coverage file that displays the non-zero coverage for each position in the genome 

```
bedtools genomecov -ibam NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram -bg > NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bedtools_cov.txt
```
Example output:    
|          |          |          |          |
|----------|----------|----------|----------|
| chr1 | 9997 | 9998  | 1 |
| chr1 | 9998 | 9999  | 2 |
| chr1 | 9999 | 10000 | 5 |


> NOTE: 
 1. column 1: chromosome or genome;
 2. column 2: start position;
 3. column 3: stop position; 
 4. column 4: read depth at the position;

#### 3. Generate a histogram of coverage data

```
bedtools genomecov -ibam NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram  > NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bedtools_hist.txt
```
Example output: 
|          |          |          |          |          |
|----------|----------|----------|----------|----------|
| chr1 | 0 | 21641911  | 248956422 | 0.0869 |
| chr1 | 1 | 7986527  | 248956422 | 0.032 |
| chr1 | 2 | 16371619 | 248956422 | 0.065 |
> NOTE: source (https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
1. column 1: chromosome (or entire genome)
2. column 2: depth of coverage from features in input file
3. column 3: number of bases on chromosome (or genome) with depth equal to column 2.
4. column 4: size of chromosome (or entire genome) in base pairs
5. column 5: fraction of bases on chromosome (or entire genome) with depth equal to column 2.

#### 4. Generate the coverage stats from the coverage file

```
awk '{sum+=$4; count++} END {print "Total number of reads across all read-covered positions: ", sum, "\nTotal number of reads-covered positions: ", count}' NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bedtools_cov.txt > NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bedtools_cov_sum_stat.txt 

awk '{ if ($4 > 0) { sum += $3 - $2 } } END { print "Total Length of Read Covered Genome Regions (bp):", sum, "\nThe proportion of GRCh38 has coverage (%):", sum*100/(3.2*10^9) }' NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bedtools_cov.txt >> NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bedtools_cov_sum_stat.txt 
```

#### 5. Coverage distribution analysis
