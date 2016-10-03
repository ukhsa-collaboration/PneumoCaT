#!/usr/bin/env python

'''
    <PneumoCaT: a tool for assigning capsular type to Streptococcus pneumoniae genomic sequence data.>
    Copyright (C) 2016 PHE

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import os, os.path, sys, subprocess, argparse, glob, yaml, inspect

module_folder_paths = ["modules"]

for module_folder_path in module_folder_paths:
  module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],module_folder_path)))
  if module_folder not in sys.path:
    sys.path.insert(1, module_folder)

import Serotype_determiner_functions
import SNP_based_Serotyping_Functions
import log_writer
from utility_functions import *


def parse_args(args):
  """
  We have set parser = argparse.ArgumentParser() and added all arguments by adding parser.add_argument.
  """
  global _parser
  
  _parser = argparse.ArgumentParser()
  _parser.add_argument('--input_directory', '-i', help='please provide the path to the directory contains the fastq files [REQUIRED - OPTION 1]')
  _parser.add_argument('--fastq_1', '-1', help='Fastq file pair 1 [REQUIRED - OPTION 2]')
  _parser.add_argument('--fastq_2', '-2', help='Fastq file pair 2 [REQUIRED - OPTION 2]')
  _parser.add_argument('--variant_database', '-d', help='variant database [OPTIONAL]; defaults to streptococcus-pneumoniae-ctvdb in the github directory', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'streptococcus-pneumoniae-ctvdb'))
  _parser.add_argument('--output_dir', '-o', help='please provide an output directory [OPTIONAL]; if none provided a pneumo_capsular_typing folder will be created in the directory containing the fastq files')
  _parser.add_argument('--bowtie', '-b', help='please provide the path for bowtie2 [OPTIONAL]; defaults to bowtie2', default='bowtie2')
  _parser.add_argument('--samtools', '-s', help='please provide the path for samtools [OPTIONAL]; defaults to samtools', default='samtools')
  _parser.add_argument('--cleanup', '-c', help='if used, all bam files generated will be removed upon completion', action = 'store_true', default='False')


  opts = _parser.parse_args(args)
  
  return opts


def main(opts):
  """
  This main function caters for three input options.

  :note option 1: if user just wants to provide the path for a dir that has the fastq files then they can with just -i option.
  :note option 2: if user chooses to provide forward and reverse files then they can specify them with -1 and -2 options.
  :note if no output dir is provided, an output dir is created in the input directory, called pneumo_capsular_typing
  """

  # This list is used to populate the fastq files later. 
  fastq_files = []
  glob_pattern = "*fastq*"
  ids = None

  # Use the utility_functions script to call check_file_exists function in common_modules which does exacly that!
  # check_file_exists(opts.bowtie2_path, 'bowtie2 path')
  # check_file_exists(opts.samtools_path, 'samtools path')

  # If an output file has not been specified, thesn create output_dir in the input_directory
  if not opts.output_dir:
    if opts.fastq_1:
      opts.output_dir =  os.path.join(os.path.dirname(opts.fastq_1), 'pneumo_capsular_typing')
      if not os.path.isdir(opts.output_dir): os.makedirs(opts.output_dir) 
    else:
      opts.output_dir = os.path.join(opts.input_directory, 'pneumo_capsular_typing')
      if not os.path.isdir(opts.output_dir): os.makedirs(opts.output_dir) #make output_directory
  else:
      if not os.path.exists(opts.output_dir):os.makedirs(opts.output_dir) 
  
  # option 1: if user wants to provide the path for a dir that has the fastq files then they can specify an input direcotry path with -i option.
  if opts.input_directory:
    check_file_exists(opts.input_directory, 'input directory')
    fastq_files = glob.glob(os.path.join(opts.input_directory, glob_pattern))

    if len(fastq_files) != 2:
      print "Unexpected number (" + str(len(fastq_files)) + ") of fastq files. Please use options -1 and -2 to specify the path to the fastq files"
      print _parser.print_help()
      sys.exit(1)

    (seqDir,seqFileName) = os.path.split(fastq_files[0])
    id,suffix = seqFileName.split(".", 1)

  # option 2: if user chooses to provide forward and reverse files then they can specify them with -1 and -2 options.
  elif opts.fastq_1 or opts.fastq_2:
    check_file_exists(opts.fastq_1, 'Fastq 1')
    check_file_exists(opts.fastq_2, 'Fastq 2')
    opts.input_directory = os.path.dirname(opts.fastq_1)
   
    fastq_files.append(opts.fastq_1)
    fastq_files.append(opts.fastq_2)
 
    (SeqDir,seqFileName) = os.path.split(fastq_files[0])  
    (id,suffix) = seqFileName.split(".",1) ## os.path.splitext(seqFileName)


  reference_fasta_file = os.path.join(opts.variant_database, 'reference.fasta')
  if opts.variant_database != 'streptococcus-pneumoniae-ctvdb':
    check_file_exists(reference_fasta_file, 'reference.fasta file')
    
  workflow = 'PneumoCaT'
  version = '1.0'
  # create a log folder in the specified input directory
  # This is set once to log all subprocesses.  The stdout and stderr log files will be in the output_dir.  The logger is called in Serotype_determiner_functions.
  if not os.path.exists(opts.output_dir + "/logs"): os.makedirs(opts.output_dir + "/logs") 
  logger = log_writer.setup_logger(info_file = opts.output_dir + "/logs/pneumo_capsular_typing.stdout", error_file = opts.output_dir + "/logs/pneumo_capsular_typing.stderr")

  #Step1: coverage based approach
  hits = Serotype_determiner_functions.find_serotype(opts.input_directory, fastq_files, reference_fasta_file, opts.output_dir, opts.bowtie, opts.samtools, opts.cleanup, id, logger, workflow=workflow, version=version) ## addition for step2
  
  ## Step2: variant based approach
  print hits
  #if len(hits) == 1 and hits[0] in ['33A', '33F']: hits = ['33A', '33F'] # force all isolates with top hit 33A or 33F to go thought the variant-based approach; 33A and 33F coverage can be right at the 90% threshold
  #elif len(hits) == 1 and hits[0] in ['11A', '11B', '11C', '11D', '11F']: hits = ['11A', '11B', '11C', '11D', '11F'] # same for the serogroup 11 isolates
  #if len(hits) > 1 or hits==['06E']:
  if hits != []:
    reference_directory = opts.variant_database
    logger = log_writer.setup_logger(info_file = opts.output_dir + "/logs/SNP_based_serotyping.stdout", error_file = opts.output_dir + "/logs/SNP_based_serotyping.stderr")
    SNP_based_Serotyping_Functions.find_serotype(opts.input_directory, fastq_files, hits, reference_directory, opts.output_dir, opts.bowtie, opts.samtools, opts.cleanup, logger, workflow=workflow, version=version)
  write_component_complete(opts.output_dir)

if __name__ == "__main__":
  opts= parse_args(sys.argv[1:])
  main(opts)
