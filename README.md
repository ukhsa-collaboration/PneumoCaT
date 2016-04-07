# README for PneumoCaT tool
---------------------------

This tool uses a two-step step approach to assign capsular type to *S.pneumoniae* genomic data. In the first step, reads from each readset are mapped to capsular locus sequences for all known capsular types (92 for S. pneumoniae plus 2 additional subtypes/molecular types). This step is considered successful if the readset matches > 90% to 1 or more capsular locus sequences. If it matches to a single capsular locus then PneumoCaT terminates here and reports this as the assigned capsular type. If more than 1 loci are matched then the tool moves to the second step; a variant based approach that utilises the capsular type variant (CTV) database to distinguish serotypes within a serogroup/genogroup. For more information you can refer to Kapatai et al 2016.

## Table of content
---------------------------

* Dependencies
* Running PneumoCaT
* PneumoCaT Output
* Examples
* Troubleshooting

## Dependencies
---------------------------

PneumoCaT  is written with Python 2.7.5 and requires the following packages installed before running:
* bowtie2/2.1.0
* samtools/0.1.19
* yaml/1.1
* numpy/python2.7/1.7.1
* lxml/python2.7.0/3.2.3
* pysam/python2.7/0.7.5
* biopython/python2.7/1.61

## Running PneumoCaT
----------------------------

```
usage: python PneumoCaT.py [-h]
                           [--input_directory INPUT_DIRECTORY]
                           [--fastq_1 FASTQ_1] [--fastq_2 FASTQ_2]

Arguments:
  -h, --help            
        show this help message and exit
  --input_directory INPUT_DIRECTORY, -i INPUT_DIRECTORY
        please provide the path to the directory contains the fastq files; [REQUIRED - OPTION 1]
  --fastq_1 FASTQ_1, -1 FASTQ_1
        Fastq file pair 1 [REQUIRED - OPTION 2]
  --fastq_2 FASTQ_2, -2 FASTQ_2
        Fastq file pair 2 [REQUIRED - OPTION 2]
  --variant_database VARIANT_DATABASE, -d VARIANT_DATABASE
        variant database [OPTIONAL]; defaults to streptococcus-pneumoniae-ctvdb in the github directory
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
        please provide an output directory [OPTIONAL]; if none provided a pneumo_capsular_typing folder will be created in the directory containing the fastq files
  --bowtie BOWTIE, -b BOWTIE
        please provide the path for bowtie2 [OPTIONAL]; defaults to bowtie2
  --samtools SAMTOOLS, -sam SAMTOOLS
        please provide the path for samtools [OPTIONAL]; defaults to samtools
```

## Output files
---------------

### STEP 1: COVERAGE-BASED APPROACH
1. **SAMPLEID.results.xml** - The XML file at step 1 is basic and it only contains the top two capsular types and their respective read coverage (% of the capsular locus length covered). If only one capsular type is matched with more than 90% coverage then the report from step 1 is considered the final result. If more than one capsular type are matched with more than 90% coverage then the software moves to step two and a second XML file is generated with the final result. Note that the output XML file from step 1 only reports two capsular types, when more could be matched and all will pass to step 2 for further distinction. If the top coverage is < 90% then no serotypes are reported and 'Failed' appears instead.
2. **SAMPLEID.sorted.bam** - BAM file generated during step 1 using the 94 capsular locus sequences as reference.
3. **SAMPLEID.sorted.bam.bai** - index file for the sorted BAM file
4. **ComponentComplete.txt** - indicates PneumoCaT analysis was completed succssfully
5. **coverage_summary.txt** - contains the coverage values for all serotypes. This is useful if the step has failed, epsecially if the top coverage falls close to the 90% threshold.

### STEP 2 - VARIANT-BASED APPROACH
1. **SAMPLEID.results.xml** - The XML file at step 2 details the assigned capsular type, total hits, the capsular types studied in this analysis as well as a full breakdown of the mutations used for this assignment. Total hits corresponds to the number of matched mutations vs the number of all mutations tested. A capsular type is assigned only if all mutations matched. If a complete match is not possible then 'Serotype Undetermined' is reported.
    Coverage statistic metrics are calculated for each gene locus used for this assignment. These include five values as detailed below:
  * Minimum depth: the minimum number of reads that maps to the gene sequence at any single position.
  * Maximum depth: the maximum number of reads that maps to the gene sequence at any single position.
  * Average depth: the average number of reads that maps across the length of the gene sequence.
  * Percent Coverage: records coverage across the length of each gene sequence.
  * Coverage Distribution: records the area of the gene sequence covered by reads.
2. **SAMPLEID.sorted.bam** - BAM file generated during step 2 using specific genes, assigned for each serogroup/genogroup in the CTV database, as reference 
3. **SAMPLEID.sorted.bam.bai** - index file for the sorted BAM file
4. **variant_summary.yml** - contains the variant matches for all serotypes tested. This is useful if a complete match was not possible and 'Serotype undetermined' was reported. It can give information on coverage for specific positions, mixed positions and other events that could lead to failure to assing a serotype.

