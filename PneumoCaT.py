#!/usr/bin/env python

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
  _parser.add_argument('--input_directory', '-i', help='please provide an input directory')
  _parser.add_argument('--workflow', '-w', help='If using a workflow you must specify an input directory')
  _parser.add_argument('--fastq_1', '-1', help='Fastq file pair 1')
  _parser.add_argument('--fastq_2', '-2', help='Fastq file pair 2')
  _parser.add_argument('--reference_fasta_file', '-se', help='serotypes directory')
  _parser.add_argument('--output_dir', '-o', help='please provide an output directory')
  _parser.add_argument('--bowtie', '-b', help='please provide the path for bowtie2', default='bowtie2')
  _parser.add_argument('--samtools', '-sam', help='please provide the path for samtools', default='samtools')
  # _parser.add_argument('--config_file_path', '-config', help='please provide the path for the config file', default='/Volumes/NGS2_DataRAID/software/workflows/config_files/determine_serotype.yml')
  opts = _parser.parse_args(args)
  
  return opts


def main(opts):
  """
  This main function caters for three input options.

  :note option 1: if user chooses to provide the workflow (e.g. streptococcus_pneumoniae.1-0), then this will look for the yaml config file and uses the path for the reference_fasta_file.
  :note option 2: if user chooses to provide forward and reverse files then they can specify them with -1 and -2 options.
  :note option 3: if user just wants to provide the path for a dir that has the fastq files then they can with just -i option.
  :note if no output dir is provided, an output dir is created in the input directory, called strep_pneumo_serotyping
  """

  # This list is used to populate the fastq files later. 
  fastq_files = []
  glob_pattern = "*.processed.*fastq*"
  ids = None
  workflow = None
  version = None

  # Use the utility_functions script to call check_file_exists function in common_modules which does exacly that!
  # check_file_exists(opts.bowtie2_path, 'bowtie2 path')
  # check_file_exists(opts.samtools_path, 'samtools path')

  # If an output file has not been specified, thesn create output_dir in the input_directory
  if not opts.output_dir:
    if opts.fastq_1:
      opts.output_dir =  os.path.join(os.path.dirname(opts.fastq_1), 'strep_pneumo_serotyping')
      if not os.path.isdir(opts.output_dir): os.makedirs(opts.output_dir) 
    else:
      opts.output_dir = os.path.join(opts.input_directory, 'strep_pneumo_serotyping')
      if not os.path.isdir(opts.output_dir): os.makedirs(opts.output_dir) #make output_directory
  else:
      if not os.path.exists(opts.output_dir):os.makedirs(opts.output_dir) 
  # This is set once to log all subprocesses.  The stdout and stderr log files will be in the output_dir.  The logger is called in Serotype_determiner_functions.
  
  # option 1: if user chooses to provide the workflow (e.g. streptococcus_pneumoniae.1-0), then this will look for the yaml config file and uses the path for the reference_fasta_file
  if opts.workflow:
    if not opts.input_directory:
      print("If using a workflow you must specify an input directory\n\n")
      print _parser.print_help()
      sys.exit(1)
    fastq_files = glob.glob(os.path.join(opts.input_directory, glob_pattern))

    # config_file = open(opts.config_file_path)
    # config = yaml.load(config_file)
    if len(fastq_files) != 2:
      print "Unexpected number (" + str(len(fastq_files)) + ") of processed fastq files"
      sys.exit(1)

    (seqDir,seqFileName) = os.path.split(fastq_files[0])
    id,suffix = seqFileName.split(".", 1)
    workflow, version = opts.workflow.split('.')
    reference_fasta_file = os.path.dirname(os.path.realpath(__file__)) + "/" + opts.workflow + "/reference.fasta"
  # option 2: if user chooses to provide forward and reverse files then they can specify them with -1 and -2 options.
  elif opts.fastq_1 or opts.fastq_2:
    check_file_exists(opts.fastq_1, 'Fastq 1')
    check_file_exists(opts.fastq_2, 'Fastq 2')
    opts.input_directory = os.path.dirname(opts.fastq_1)
    if not opts.reference_fasta_file:
      print("If you are using -1 and -2 options, you must provide the serotype directory folder")
      print _parser.print_help()
      sys.exit(1)
    # not using line 103 here as no input_directory is provided when supplied with -1 and -2 so it doesnt know where to create the output directory.
    
    fastq_files.append(opts.fastq_1)
    fastq_files.append(opts.fastq_2)
    reference_fasta_file = os.path.join(opts.reference_fasta_file, 'reference.fasta')
    workflow, version = os.path.basename(opts.reference_fasta_file).split('.')
    (SeqDir,seqFileName) = os.path.split(fastq_files[0])  
    (id,ext) = seqFileName.split(".",1) ## os.path.splitext(seqFileName)

  # option 3: if user just wants to provide the path for a dir that has the fastq files then they can with just -i option.
  elif opts.input_directory:
    check_file_exists(opts.input_directory, 'input directory')
    if not opts.reference_fasta_file:
      print("If you are providing the input directory for the fastqs, you must provide the serotype directory folder")
      print _parser.print_help()
      sys.exit(1)
    fastq_files = glob.glob(opts.input_directory + "/" + glob_pattern)
    if fastq_files == []: fastq_files = glob.glob(os.path.join(opts.input_directory, '*.fastq*'))
    reference_fasta_file = os.path.join(opts.reference_fasta_file, 'reference.fasta')
    workflow, version = os.path.basename(opts.reference_fasta_file).split('.')
    (seqDir,seqFileName) = os.path.split(fastq_files[0])  
    (id,ext) = seqFileName.split(".",1) ## os.path.splitext(seqFileName)
  if not os.path.exists(opts.input_directory + "/logs"): os.makedirs(opts.input_directory + "/logs") # if script is run outside of the pipeline, the logs folder is not present and the script will exit unless the logs folder is created... 
  logger = log_writer.setup_logger(info_file = opts.input_directory + "/logs/strep_pneumo_serotyping.stdout", error_file = opts.input_directory + "/logs/strep_pneumo_serotyping.stderr")

  hits = Serotype_determiner_functions.find_serotype(opts.input_directory, fastq_files, reference_fasta_file, opts.output_dir, opts.bowtie, opts.samtools, id, logger, workflow=workflow, version=version) ## addition for step2
  ## addition for step2: SNP_based_serotyping
  print hits
  if len(hits) == 1 and hits[0] in ['33A', '33F']: hits = ['33A', '33F']
  elif len(hits) == 1 and hits[0] in ['11A', '11B', '11C', '11D', '11F']: hits = ['11A', '11B', '11C', '11D', '11F']
  if len(hits) > 1 or hits==['06E']:
    reference_directory = os.path.dirname(reference_fasta_file)
    logger = log_writer.setup_logger(info_file = opts.input_directory + "/logs/SNP_based_serotyping.stdout", error_file = opts.input_directory + "/logs/SNP_based_serotyping.stderr")
    SNP_based_Serotyping_Functions.find_serotype(opts.input_directory, hits, reference_directory, opts.output_dir, opts.bowtie, opts.samtools, logger, workflow=workflow, version=version)
  write_component_complete(opts.output_dir)

if __name__ == "__main__":
  opts= parse_args(sys.argv[1:])
  main(opts)
