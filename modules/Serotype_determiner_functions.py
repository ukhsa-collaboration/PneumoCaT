"""

.. module:: Serotype_determiner_functions.py

.. moduleauthor:: Ali Al-Shahib, Georgia Kapatai

"""
import os, os.path, sys, subprocess, inspect

module_folder_paths = ["modules"]

for module_folder_path in module_folder_paths:
  module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],module_folder_path)))
  if module_folder not in sys.path:
    sys.path.insert(1, module_folder)
import pysam
import glob
import string, re
import fileinput
import shutil
import numpy
from Bio import SeqIO
from lxml import etree
import log_writer
from utility_functions import *
import fnmatch
from os.path import exists
import time
import yaml, operator



def find_serotype(input_directory, fastqs, reference_fasta_file, output_dir, bowtie, samtools, clean_bam, id, logger, workflow = "", version = ""):
  
  """
  
  This is the main function that calls other functions to find serotype.

  :param fastqs: directory that contains two fastq files.  It will be in the following format: id.workflow.version.suffix, e.g. 1.strep_pneumo.1_1.1.trimmed.fastq
  :type fastqs: directory
  :param reference_fasta_file: reference_fasta_fileerence file to be used in the mapping function.  In this case it should be a fasta file that contains the serotypes.
  :param output_dir: Output directory 
  :type output_dir: directory 
  :param bowtie: path to bowtie
  :type bowtie: path
  :param samtools: path to samtools
  :type samtools: path

  :returns: Top two serotypes with their percentage coverage.
  :rtype: string
  :rtype: float

  """
  bam = try_and_except(output_dir + "/logs/strep_pneumo_serotyping.stderr", mapping, input_directory, fastqs, reference_fasta_file, output_dir, bowtie, samtools, id, logger)
  output_file = open(output_dir + "/" + id + ".results.xml", "w")
  try_and_except(output_dir + "/logs/strep_pneumo_serotyping.stderr", best_coverage,bam,reference_fasta_file)
  hits = try_and_except(output_dir + "/logs/strep_pneumo_serotyping.stderr", output_all,bam,reference_fasta_file,output_file,id,workflow,version) # added for step 2
  try_and_except(output_dir + "/logs/strep_pneumo_serotyping.stderr", cleanup,output_dir)
  if clean_bam == True:
    bamfiles = glob.glob(bam+'*')
    for file in bamfiles: os.remove(file)
  return hits

def mapping(input_directory, fastqs, reference_fasta_file_path, output_dir, bowtie, samtools, id, logger):

  """

  This function runs bowtie for mapping of the fastq files with the reference_fasta_fileerence file provided.

  :param fastqs: directory that contains two fastq files.  It will be in the following format: id.workflow.version.suffix, e.g. 1.strep_pneumo.1_1.1.trimmed.fastq
  :type fastqs: directory
  :param reference_fasta_file: reference_fasta_file, this is defined in the find_serotype function.
  :type reference_fasta_file: file
  :param bowtie: path to bowtie
  :type bowtie: path
  :param samtools: path to samtools
  :type samtools: path

  :returns: sorted and indexed bam file.
  :rtype: file

  """
  try:
    os.makedirs(output_dir + "/tmp")
  except OSError:
    if os.path.isdir(output_dir + "/tmp"):
      # We are nearly safe
      pass
    else:
      # There was an error on creation, so make sure we know about it
      raise

  null = open(os.devnull, 'w')

  bam_sorted = os.path.join(output_dir, id + '-sorted')
  bam_out = os.path.join(output_dir, id + '-sorted.bam')
  sam_parsed = os.path.join(output_dir + "/tmp", id + '.tmp') # temporary sam output
  sam = os.path.join(output_dir + "/tmp", id + '.sam')
  bam = os.path.join(output_dir + "/tmp", id + '.bam')
  
 
  # copy the reference fasta file to the tmp directory and index
  reference_fasta_file = output_dir + "/tmp/reference.fasta"
  shutil.copyfile(reference_fasta_file_path, reference_fasta_file)
  print "running bowtie index"
  bowtie_index=  bowtie + "-build"
  log_writer.info_header(logger, "Creating reference_fasta_fileerence index")
  process = subprocess.Popen([bowtie_index, reference_fasta_file, reference_fasta_file], stderr=subprocess.PIPE, stdout=subprocess.PIPE) # generate index of reference_fasta_fileerence fasta for mapping
  process.wait()
  log_writer.log_process(logger, process)
  
  # # run bowtie
  cmd = [bowtie]
  cmd += [ '--fr', '--minins', '300', '--maxins', '1100', '-x', reference_fasta_file, '-1', fastqs[0], '-2', fastqs[1],'-S', sam, '-k', '99999', '-D', '20', '-R', '3', '-N', '0', '-L', '20', '-i', 'S,1,0.50'] # write to tmp 
  log_writer.info_header(logger, "Running bowtie to generate sam file")
#  print "running bowtie"
  process = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
  process.wait()
  log_writer.log_process(logger, process, log_error_to = "info")

  # run remove_secondary_mapping_bit to deduct 256 from any bit that is above 256 using the sam file.  The sam_parsed file is the output that is used to convert to bam.

  try_and_except(input_directory + "/logs/strep_pneumo_serotyping.stderr", remove_secondary_mapping_bit, sam, sam_parsed)

  log_writer.info_header(logger, "Convert sam to bam")
  process = subprocess.Popen([samtools, 'view', '-bhS', '-o', bam, sam_parsed], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
  process.wait()
  log_writer.log_process(logger, process, log_error_to = "info")
  # sort bam
  log_writer.info_header(logger, "Sort the bam file")
  process = subprocess.Popen([samtools, 'sort', bam, bam_sorted], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
  process.wait()
  log_writer.log_process(logger, process, log_error_to = "info")
  # index bam
  log_writer.info_header(logger, "Index the BAM file")
  process = subprocess.Popen([samtools, 'index', bam_sorted + ".bam"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
  process.wait()
  log_writer.log_process(logger, process, log_error_to = "info")
  return bam_out 

def remove_secondary_mapping_bit(sam,sam_parsed):

  """
  
  This function removes secondary mapping bit from the FLAG column in SAM by deducting 256 from any bit that is above 256.

  :param sam: sam file that has been produced from the mapping function
  :type sam: file
  :param sam_parsed: the new parsed sam file.
  :type sam_parsed: file

  :returns: parsed sam file ready to be converted to bam.
  :rtype: file

  """
  
  lines = iter(fileinput.input([sam]))
  sam_parsed_file = open(sam_parsed, "w")
  for line in lines:
    if line.startswith('@'):
      sam_parsed_file.write(line)
    else:
      # chomp line
      line = line.rstrip('\n')
      details = line.split("\t")
      flag = int(details[1])
      if flag > 256:
        details[1] = str(flag - 256)
      print >> sam_parsed_file, '\t'.join(details)
        # print headers
  sam_parsed_file.close()

def reference_fasta_fileerence_length(reference_fasta_file):

  """
  The feature_length function finds the length and id of each contig in a fasta file.

  :param reference_fasta_file: The reference_fasta_fileerence fasta file (if used for finding serotypes, then each contig in the fasta file should be a serotype)
  :type reference_fasta_file: file

  :returns: lengths  - list of total lengths of each serotype
  :rtype: list
  :returns: features_id -  list of serotype ids
  :rtype: list

  """

  lengths = []
  ids = []
  for fasta_record in SeqIO.parse(open(reference_fasta_file,"r"), "fasta"):
    lengths.append(len(fasta_record.seq))
    ids.append(fasta_record.id)
  return (lengths,ids)


def best_coverage(bam,fasta):

  """
  The best_coverage function finds the serotype that has the highest read coverage over the contig (or serotype) length and assigns that as the serotype for the particular strain.

  :param bam: BAM file of the strain
  :type bam: file
  :param fasta: Fasta file of serotypes
  :type fasta: file

  :returns: serotype - The correct serotype for the strain.
  :rtype: string
  :returns: percentage coverage - percentage of reads covering the length of the top two serotypes
  :rtype: float

  """

  # NOTE: id HERE IS SEROTYPE id (e.g. 06A) NOT NGS id
  # use pysam to read bam file
  bamfile = pysam.Samfile(bam, "rb")
  highest_percent_coverage = 0
  second_highest_coverage = 0
  selected_id = ""
  second_selected_id = ""
  coverage = {} ## addition for step2
  mean_depth = {}
  # get lengths and ids from reference_fasta_fileerence_length function
  (lengths,ids) = reference_fasta_fileerence_length(fasta)
  for determine_reference_fasta_fileerence_length,id in zip(lengths,ids):
    bases_with_coverage_array = []
    all_bases_coverage = []
    # go through each id and perform pileup
    for pileupcolumn in bamfile.pileup(id):
      all_bases_coverage.append(pileupcolumn.n)
    # Any coverage greater than zero
      if pileupcolumn.n > 4:
        # append all base coverages for the id in a list 
        bases_with_coverage_array.append(pileupcolumn.n)
    # find the overall percentage coverage for each id
    percent_coverage_for_id = len(bases_with_coverage_array)/float(determine_reference_fasta_fileerence_length+1)*100
    mean_depth[id] = sum(all_bases_coverage)/float(len(all_bases_coverage)) if len(all_bases_coverage)>0 else 0
    # As the highest percentage coverage is set to zero at the start of the script, this will increase as the script loops through all ids (or serotypes). Until it finds no percentage higher than the previous one, it assigns the highest_percent_coverage = percent_coverage_for_id.  Same story for the second highest coverage, its the one before the highest.
    coverage[id] = percent_coverage_for_id ## addition for step2
    if percent_coverage_for_id > highest_percent_coverage:
      second_highest_coverage = highest_percent_coverage
      second_selected_id = selected_id
      highest_percent_coverage = percent_coverage_for_id
      selected_id = id 
    elif percent_coverage_for_id > second_highest_coverage:
      second_highest_coverage = percent_coverage_for_id
      second_selected_id = id
  
  ## export hits into a yaml file to investigate failed samples
  out_fp = open(os.path.join(os.path.dirname(bam), 'coverage_summary.txt'), 'w')
  out_fp.write('Serotype\tCoverage\tDepth\n')
  sorted_coverage = sorted(coverage.items(), key=operator.itemgetter(1), reverse=True)
  for ser, cov in sorted_coverage:
    out_fp.write(ser+'\t'+str(cov)+'\t'+str(mean_depth[ser])+'\n')
  out_fp.close()

  ## select hits > 90%
  hits = [f for f in coverage.keys() if coverage[f]>90] ## addition for step2
  # return these so it is used by the output_all function
  return (selected_id, second_selected_id, highest_percent_coverage, second_highest_coverage, hits) ## addition for step2

  bamfile.close


def output_all(bam,fasta,output_file,ngs_sample_id,workflow="",version=""):

  """
  This function calculates the min_depth, meanQ and mean_depth and outputs all the required results in xml format.

  :param bam: BAM file of the strain
  :type bam: file
  :param fasta: Fasta file of serotypes
  :type fasta: file
  :param output_file
  :type output_file: file
  :param ngs_sample_id: Taken from mapping function
  :type ngs_sample_id: int
  :param workflow: Taken from mapping function
  :type workflow: string
  :param version: Taken from mapping function
  :type version: int

  :returns: serotype - The correct serotype for the strain.
  :rtype: string
  :returns: highest percentage coverage - percentage of reads covering the length of the top serotype
  :rtype: float
  :returns: QC_mean_depth - mean depth of the highest percentage coverage
  :rtype: int
  :returns: QC_minimum_depth - min depth of the highest percentage coverage
  :rtype: int
  :return: QC_meanQ - mean base quality of the highest percentage coverage 
  :rtype: int
  :returns: second_serotype - The second most probable serotype for the strain.
  :rtype: string
  :return: second_highest_coverage - percentage of reads covering the length of the second most probable serotype for the strain.
  rtype: float


  """

  (selected_id, second_selected_id, highest_percent_coverage, second_highest_coverage, hits) = best_coverage(bam,fasta) ## addition for step2

  bamfile = pysam.Samfile(bam, "rb")
  quals = []
  overall_qual = []
  bases_with_coverage_array = []
  for pileupcolumn in bamfile.pileup(selected_id):
    for pileupread in pileupcolumn.pileups:
      # qual_ascii is the quality for each base but in ascii format.
      qual_ascii = pileupread.alignment.qual[pileupread.qpos]
      # phred score is an integer representing the Unicode code point of the character in qual_ascii
      phred_score = ord(qual_ascii) - 33 
      # append all quals in one whole list
      quals.append(phred_score)
      if pileupcolumn.n > 4:
        # append all base coverages for the id in a list 
        # This is done to find max and mean depth
        bases_with_coverage_array.append(pileupcolumn.n)

  # mean base quality
  
  
  if not bases_with_coverage_array:
    selected_id = "Failed"
    second_selected_id = "Failed"
    min_depth = 0
    meanQ = 0
    mean_depth = 0
  else:
    min_depth = min(bases_with_coverage_array)
    meanQ = round(numpy.mean(quals),1)
    mean_depth = round(numpy.mean(bases_with_coverage_array),1)

  if highest_percent_coverage <= 90 or mean_depth <= 20:
    selected_id = "Failed"
    second_selected_id = "Failed"
  
  root = etree.Element("ngs_sample", id = ngs_sample_id)
  # if a workflow is provided, i.e. for the service, then the xml output will contain a workflow.  If not (e.g. when you use -1 and -2 options) then it will skip it as workflow is not provided.
  if workflow:
    workflow_xml = etree.SubElement(root, "workflow", value= workflow, version=version)

  results = etree.SubElement(root, 'results')
  result = etree.SubElement(results, "result", type="Serotype", value = selected_id)
  etree.SubElement(result, "result_data", type='QC_coverage', value = str("%.2f" % highest_percent_coverage)) # only two decimal places.Also float -> str
  etree.SubElement(result,"result_data", type='QC_mean_depth', value =  str(mean_depth))
  etree.SubElement(result, "result_data", type='QC_minimum_depth', value = str(min_depth))
  etree.SubElement(result, "result_data", type='QC_meanQ', value = str(meanQ))
  etree.SubElement(result, "result_data", type='second_value', value = second_selected_id)
  etree.SubElement(result, "result_data", type='QC_coverage_second_value', value = str("%.2f" % second_highest_coverage))
#  print etree.tostring(root, pretty_print=True)

  print >> output_file, etree.tostring(root, pretty_print=True)

  bamfile.close()
  output_file.close()
  ## addition for step2
  if selected_id != 'Failed':
    return hits
  else:
    return [] 

def cleanup(output_dir):

  # remove tmp file in output_dir

  try:
    shutil.rmtree(output_dir + "/tmp/")
  except OSError:
    raise
