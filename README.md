# README for PneumoCaT tool
---------------------------

This tool uses a two-step step approach to assign capsular type to S.pneumoniae genomic data. In the first step, reads from each readset are mapped to capsular locus sequences for all known capsular types (92 for S. pneumoniae plus 2 additional subtypes/molecular types). This step is considered successful if the readset matches > 90% to 1 or more capsular locus sequences. If it matches to a single capsular locus then PneumoCaT terminates here and reports this as the assigned capsular type. If more than 1 loci are matched then the tool moves to the second step; a variant based approach that utilises the capsular type variant (CTV) database to distinguish serotypes within a serogroup/genogroup. For more information you can refer to Kapatai et al 2016.

## Table of content
---------------------------

* Dependencies
* Running PneumoCaT
* PneumoCaT Output
* Examples

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
