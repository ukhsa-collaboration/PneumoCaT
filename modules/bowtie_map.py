#!/usr/bin/env python

"""
.. module:: bowtie_map.py
  :platform: Unix, Windows

.. moduleauthor:: Georgia Kapatai

"""

import os, subprocess, re, glob, inspect, sys, argparse
import datetime
import shutil

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

module_folder_paths = ["modules", '~/common_modules']
for module_folder_path in module_folder_paths:
    module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe()))[0],module_folder_path)))
    if module_folder not in sys.path:
        sys.path.insert(1, module_folder)
import log_writer
from utility_functions import *

# allow running of main automatically
#from automain import *


def create_index(reference_fasta_file, bowtie_path, samtools_path, tmpdir, logger):
    'Build a bowtie2 index and fai index from the given input(s)'

    fasta_file = tmpdir + "reference.fasta"
    shutil.copyfile(reference_fasta_file, fasta_file)
    bt2_index = fasta_file + '.1.bt2'
    fai_index = fasta_file + '.fai'
    log_writer.info_header(logger, 'Building bowtie2 index for {}...'.format(fasta_file))
    if os.path.exists(bt2_index):
        log_writer.write_log(logger, 'Bowtie2 index for {} is already built...'.format(fasta_file), 'info')
    else:
        bowtie2_index = bowtie_path+'-build'
        subprocess.call([bowtie2_index, '-f', fasta_file, fasta_file])
    
    log_writer.info_header(logger, 'Building samtools index for {}...'.format(fasta_file))
    if os.path.exists(fai_index):
        log_writer.write_log(logger, 'Samtools index for {} is already built...'.format(fasta_file), 'info')
    else:
        subprocess.call([samtools_path, 'faidx', fasta_file])

    return fasta_file

def modify_bowtie_sam(samfile, logger):
    'Modify SAM formatted output from Bowtie to maintain secondary alignments for downstream pileup'
    with open(samfile) as sam, open(samfile + '.mod', 'w') as sam_mod:
        
        log_writer.info_header(logger, 'Modifying SAM formatted output from Bowtie to maintain secondary alignments for downstream pileup...')
        for line in sam:
            if not line.startswith('@'):
                fields = line.split('\t')
                flag = int(fields[1])
                flag = (flag - 256) if (flag > 256) else flag
                sam_mod.write('\t'.join([fields[0], str(flag)] + fields[2:]))
            else:
                sam_mod.write(line)
    return samfile + '.mod'

def clean_up(tmp, logger):
    ' Remove temporary files'
    log_writer.info_header(logger, 'Removing temporary files:'+tmp)
    try:
        shutil.rmtree(tmp)
    except OSError:
        raise


def mapping(fastq_file, reference, bowtie_path, samtools_path, outdir, logger):
    prefix =  re.sub('.R\d{1}.processed.fastq.?gz?|.processed.R\d{1}.fastq.?gz?', '', os.path.basename(fastq_file))
    sample = os.path.basename(fastq_file).split('.')[0]
    # create index if not available
    tmp = outdir + "/{0}_tmp/".format(sample)
    if not os.path.exists(tmp): os.mkdir(tmp)
    reference_fasta_file = create_index(reference, bowtie_path, samtools_path, tmp, logger)

    fastq1 = fastq_file
    #if _args.source == 0:
    fastq2 = fastq_file.replace('.R1', '.R2') if fastq_file.find('.R1') != -1 else fastq_file.replace('_1.processed', '_2.processed')
    #else:
        #fastq2 = fastq_file.replace('_1.fq', '_2.fq')
    samfile =  tmp + prefix + '.sam'
    bamfile = tmp + prefix + '.bam'
    sorted_bam_prefix = outdir + '/' + prefix + '.sorted'
    if not os.path.isfile(sorted_bam_prefix+'.bam'): # change this to bamfiles/prefix.sorted.bam
        # create sam file
        log_writer.info_header(logger, "Running bowtie to generate sam file")
        subprocess.call([bowtie_path, '--fr', '--minins', '300', '--maxins', '1100', '-x', reference_fasta_file, '-1', fastq1, '-2', fastq2, '-S', samfile, '-k', '99999', '-D', '20', '-R', '3', '-N', '0', '-L', '20', '-i', 'S,1,0.50']) # write to tmp

        # remove flags > 256 to allow reads to map in more than one locations
        sam_mod = modify_bowtie_sam(samfile, logger)

        # convert to bam
        log_writer.info_header(logger, "Running samtools to convert sam to bam file")
        subprocess.call([samtools_path, 'view', '-bS', '-o', bamfile, sam_mod])

        # sort bam file
        log_writer.info_header(logger, "Sort the bam file")
        subprocess.call([samtools_path, 'sort', bamfile, sorted_bam_prefix])

        # index bam file
        log_writer.info_header(logger, "Index the BAM file")
        subprocess.call([samtools_path, 'index', sorted_bam_prefix+'.bam'])

        # get stats
        output = subprocess.check_output([samtools_path,  "flagstat", sorted_bam_prefix+'.bam' ])

        ## clean up all unnecessary files
        #clean_up(tmp, logger)

    return sorted_bam_prefix+'.bam', reference_fasta_file

def pileup(sorted_bamfile, reference, samtools, outdir, logger):
    ' Create a mpileup file '
    log_writer.info_header(logger, 'Create mpileup file')
    filename = os.path.join(outdir, os.path.basename(sorted_bamfile).split('.')[0] + '.mpileup')
    pileupFile = open(filename, 'w')
    process = subprocess.Popen([samtools, 'mpileup', '-B', '-A', '-f', reference, sorted_bamfile], stderr=subprocess.PIPE, stdout=subprocess.PIPE) # -A -count anomalous read pairs, -B - disable BAQ computation  and -f FILE - indexed reference sequence file
    result = process.stdout
    for l in result:
        pileupFile.write(l)
    pileupFile.close()


