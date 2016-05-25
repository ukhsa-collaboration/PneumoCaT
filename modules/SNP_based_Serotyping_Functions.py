#!/usr/bin/env python

import os, sys, subprocess, argparse, glob, yaml, inspect, re
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import pickle, datetime
import shutil
from lxml import etree

module_folder_paths = ["modules", '~/common_modules']
for module_folder_path in module_folder_paths:
    module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe()))[0],module_folder_path)))
    if module_folder not in sys.path:
        sys.path.insert(1, module_folder)

import bowtie_map
from operator import itemgetter
from itertools import groupby
from collections import Counter
from Bio import pairwise2
from utility_functions import *
import log_writer
import yaml

# ============================================================= #
# Small functions                                               #
# ============================================================= #

def clean_up(tmp, logger):
    ' Remove temporary files'
    log_writer.info_header(logger, 'Removing temporary files:'+tmp)
    try:
        shutil.rmtree(tmp)
    except OSError:
        raise
    
def parse_reference_fasta(fasta_file):
    """
    Function to parse the reference fasta file provided by function parse_serotype_specific_info.
    
    Arguments:
    fasta_file: path to reference fasta file
    """
    refseq = {}
    for record in SeqIO.parse(open(fasta_file, 'rU'), 'fasta'):
        refseq[record.id] = record.seq

    return refseq


def parse_serotype_specific_info(serotypes, reference_directory, output_dir, fastq, logger, workflow="", version=""):
    """
    Function that provides the path for reference fasta file and the mutation database for the
    group of serotypes provided. The reference information for all genogroups is stored in folder
    provided in this function as data_path variable. Within this folder each genogroup has a
    subfolder that contains the reference.fasta file and the mutation database saved in .pickled
    and .yaml files.

    Arguments:
    serotypes: a list provided by Serotype_determiner_functions.find_serotype() function
    output_dir: the path to the output directory
    fastq: path to the fastq file.
    """
    global _serotypes
    # hard code the data_path variable to the reference folder.
    data_path = reference_directory
    _serotypes = sorted(serotypes)
    key = '_'.join(_serotypes)
    datadir = os.path.join(data_path, key)
    if not os.path.exists(datadir):
        # if the path to the data directory is not found then there are two posibilities:
        # a) not all serotypes passed the coverage > 90% thershold of Serotype_determiner_functions.find_serotype() function
        # b) the group of serotypes provided correspond to a combination that has not been encounter before; either because
        # it's a mixed culture or a genogroup that has no data available i.e. serogroups 24 or 32.
        # The following line interrogates all data directories to check whether the serotypes provided
        # correspond to a partial genogroup (option a) or not (option b)
        k = [f for f in os.listdir(data_path) if set(_serotypes)-set(f.split('_')) == set([])]
        if k:
            key = k[0]
            datadir = os.path.join(data_path, key)
        else:
            if len(_serotypes) == 1:
                print 'Analysis completed'
                clean_up(output_dir, logger)
            else:   
                print "No reference data available", str(serotypes)
                special_output_xml(serotypes,output_dir, fastq, logger, workflow, version)
            #sys.exit(1)
            unique = None
            reference_fasta_file = None
            return unique, reference_fasta_file, _serotypes
    data_file = open(os.path.join(datadir,'mutationdb.pickled'), 'rb')
    unique = pickle.load(data_file)
    reference_fasta_file = os.path.join(datadir, 'reference.fasta')
    serotypes = os.path.basename(datadir).split('_')

    return unique, reference_fasta_file, serotypes

def special_output_xml(serotypes, output_dir, filepath, logger, workflow="", version=""):
    """
    Function for writing a results.xml file if no reference data are available for the serotypes
    provided. This could be either due to a mixed sample or if the sample belongs to the two
    serogroups (24 and 32) that have no reference data available yet.

    Arguments:
    serotypes: list with the serotypes provided in the wrapper script strep_pneumo_serotyping.py
    output_dir: path tp the output directory
    filepath: path to the fastq R1 file
    workflow: pipeline workflow provided in the wrapper script strep_pneumo_serotyping.py; optional argument
    version: version for pipeline workflow provided in the wrapper script strep_pneumo_serotyping.py; optional argument
    """ 
    ngs_sample_id = os.path.basename(filepath).split('.',1)[0]
    serotypes.sort()
    output_file = open(output_dir + "/" + ngs_sample_id + ".results.xml", "w")

    if serotypes in [['24B'], ['24F'], ['24B','24F'], ['24A','24B','24F']]:
        result = 'Reference data not available: Serogroup 24'
    elif serotypes == ['32A','32F']:
        result = "Reference data not available: Serogroup 32"
    else:
        result = 'Mixed'+str(serotypes)

    root = etree.Element("ngs_sample", id = ngs_sample_id)
    if workflow:
        workflow_xml = etree.SubElement(root, "workflow", value= workflow, version=version)
    results = etree.SubElement(root, 'results')
    comment = etree.Comment('(START) Serotype Distinction Results (START)')
    results.append(comment)
    result = etree.SubElement(results, "result", type="Serotype_Distinction", value = result)
    comment = etree.Comment('(END) Serotype Distinction Results (END)')
    results.append(comment)
    print >> output_file, etree.tostring(root, pretty_print=True)
    output_file.close()

def output_xml(match, output_dir, ngs_sample_id, unique, logger, workflow="", version=""):
    """
    Function for writing the results.xml file reporting the matched serotype and the mutations that were used to
    predict said serotype.

    Arguments:
    match: list generated by the function pileup_investigation. Contains a list [serotype, mcnt, tcnt, mutations]
           for each serotype
    output_dir: path tp the output directory
    ngs_sample_id: sample name
    unique: dictionary with all mutations tested for this serogroup
    workflow: pipeline workflow provided in the wrapper script strep_pneumo_serotyping.py; optional argument
    version: version for pipeline workflow provided in the wrapper script strep_pneumo_serotyping.py; optional argument
    """
    output_file = open(output_dir + "/" + ngs_sample_id + ".results.xml", "w")
    serotypes = [m[0] for m in match]
    if [f for f in serotypes if f.startswith('06E')]:
        serotypes = [f for f in serotypes if not f.startswith('06E')] + ['06E']
    serotypes = sorted(serotypes)
    serotypes_tested = ','.join(serotypes)
    scores = [(f[0], f[1]/float(f[2])) for f in match] #creates a list of tuples [(serotype1, score1), (serotype2, score2), ...] where score = matched_hits/total_mutations_tested
    scores = sorted(scores, key=lambda x: x[1], reverse=True) # sorts larger to smaller based on scores
    top_hit, top_score = scores[0]#
    mixed_serotypes = ''
    for record in match:
        if record[0] == top_hit:
            serotype, mcnt, tcnt, mutations = record
            if '06E' in serotype: serotype='06E'
    if mcnt != tcnt and (serotypes != ['15A','15B', '15C', '15F'] and serotype != '06E'): # reports a failure if the top hit scored < 1.
        if [f for f in mutations if f[3].startswith('Mixed:')]:
            if len([f for f in scores if f[1] == top_score]) > 1:
                top_hit = 'Mixed: ' + str([f[0] for f in scores if f[1] == top_score])
            else:
                top_hit = top_hit + '+'
        else:
            if len(_serotypes) == 1 and top_hit == _serotypes[0]:
                top_hit = _serotypes[0]+'+'
            else:
                top_hit = 'Serotype undetermined'
                
        # even though no serotype can be determined failure at this stage will still report the mutations for the top hit.
    elif mcnt != tcnt and serotypes == ['15A','15B', '15C', '15F']:
        if [f for f in mutations if f[0] =='wciZ' and f[3].startswith('Mixed:')]: # check if mixed for wciZ
            mcnt += 1
            if mcnt == tcnt:
                top_hit = '15B/C'
            else:
                if [f for f in mutations if f[0] !='wciZ' and f[3].startswith('Mixed:')] and len([f for f in scores if f[1] == top_score]) > 1:
                    top_hit = 'Mixed: ' + str(['15A','15B/C'])
                else:
                    top_hit = 'Serotype undetermined'
        else:
            if len(_serotypes) > 1 and [f for f in mutations if f[3].startswith('Mixed:')]:
                top_hit = top_hit+'+' if len([f for f in scores if f[1] == top_score]) == 1 else 'Mixed: '+ str([f[0] for f in scores if f[1] == top_score])
            else:
                if len(_serotypes) == 1 and top_hit == _serotypes[0]:
                    top_hit = _serotypes[0]+'+'
                else:
                    top_hit = 'Serotype undetermined'
    elif len([f[1] for f in scores if f[1] == 1])> 1: # reports mixed sample if more than one serotypes scored == 1
        mixed_serotypes = [f[0] for f in scores if f[1] == 1]
        top_hit = 'Mixed: ' + str([f[0] for f in scores if f[1] == 1])
    #elif top_score<1 and len([f[1] for f in scores if f[1] == top_score]) > 1:
    #    mixed_serotypes = [f[0] for f in scores if f[1] == top_score]
    #    top_hit = 'Mixed: ' + str(mixed_serotypes) 
    ## create xml ##
    root = etree.Element("ngs_sample", id = ngs_sample_id)
    # if a workflow is provided, i.e. for the service, then the xml output will contain a workflow.If not (e.g. when you use -1 and -2 options) then it will skip it as workflow is not provided.
    if workflow: # if a workflow is provided, i.e. for the service, then the xml output will contain a workflow otherwise no workflow will be reported.
        workflow_xml = etree.SubElement(root, "workflow", value= workflow, version=version)

    results = etree.SubElement(root, 'results')
    comment = etree.Comment('(START) Serotype Distinction Results (START)')
    results.append(comment)
    result = etree.SubElement(results, "result", type="Serotype_Distinction", value = str(top_hit))
    if mixed_serotypes:
        comment = etree.Comment('(END) Serotype Distinction Results (END)')
        results.append(comment)
        print >> output_file, etree.tostring(root, pretty_print=True)
        output_file.close()
    else:
        result = etree.SubElement(results, "result", type="Serotype_Distinction_Total_Hits", value = '{0}/{1}'.format(mcnt, tcnt))
        result = etree.SubElement(results, "result", type="Serotype_Distinction_Serotypes_Tested", value = serotypes_tested)
        for key in unique.keys():
            genes = unique[key].keys() if serotype not in unique[key].keys() else unique[key][serotype].keys()
            for gene in genes:
                # extract necessary information
                mutcnt = len(unique[key][gene]) if key == 'snps' else 1
                gene_specific_mutations = [m for m in mutations if m[0] == gene and m[4] == key]
                if gene_specific_mutations == []: continue
                result = etree.SubElement(results, "result", type="Serotype_Distinction_Gene", value = gene)
                hits = '{0}/{1}'.format(len([m for m in gene_specific_mutations if m[1] != 'Failed']),mutcnt)
                failed = [m for m in gene_specific_mutations if m[1] == 'Failed']
                if key in ['genes', 'pseudo']:
                    coverage, avg_depth, min_depth, max_depth, ranges, meanQ, length = _metrics[gene]
                elif key == 'snps':
                    if serotype == '06E' and gene == 'wciP':
                        st = top_hit.split('(')[1].split(')')[0]
                        coverage, avg_depth, min_depth, max_depth, ranges, meanQ, length = _metrics[gene][st]
                    else:
                        coverage, avg_depth, min_depth, max_depth, ranges, meanQ, length = _metrics[gene][serotype]
                else:
                    coverage, avg_depth, min_depth, max_depth, ranges, meanQ, length = _metrics[gene][unique['allele'][gene][serotype]]
                cov_distribution = ";".join(['{0}-{1}'.format(f[0], f[1]) for f in ranges])
                failures = [m[3] for m in gene_specific_mutations if m[3]!='None']
                failure_tag = ";".join(failures) if failures!=[] else "None"
                detected = 'N' if [f for f in gene_specific_mutations if (f[1] == 'Failed' and coverage < 0.8 and key!='genes') or (key == 'genes' and f[1] == 0)] else 'Y'
                result_tags = {}
                key_tag = {'SNPs':'snps', 'Pseudogene':'pseudo', 'Allele':'allele'}
                # Gather results requires different approach based on the mutation tag for each gene 
                if key == 'snps': # reports snps using format 'nt_pos codon: aa_pos aa'
                    snps = ['{0} {1}:{2} {3}'.format(f[1], f[2][0], f[2][1], f[2][2]) for f in gene_specific_mutations if f[2] != 'Failed']
                    result_tag = ";".join(snps) if snps!=[] else "Failed"
                    if failed == []: hits = '{0}/{1}'.format(len(snps),mutcnt)
                elif key == 'allele': # reports allele name 
                    result_tag = gene_specific_mutations[0][1]
                elif key == 'pseudo':
                    # based on serotype low coverage could be consider failure if gene is supposed to be functional in tested serotype
                    if gene_specific_mutations[0][1] == 'Failed':
                        result_tag = 'Failed'
                    else:
                        # or a pseudogene if this is expected for tested serotype
                        result_tag = 'N' if gene_specific_mutations[0][1]==0 else gene_specific_mutations[0][2].replace('Pseudogene','Y')
                etree.SubElement(result, "result_data", type='Detected', value = detected) 
                etree.SubElement(result, "result_data", type='Hits', value =hits)
                for tag in key_tag.keys():
                    if key_tag[tag] == key:
                        etree.SubElement(result, "result_data", type=tag, value=result_tag)
                    else:
                        etree.SubElement(result, "result_data", type=tag, value="")
                etree.SubElement(result, "result_data", type='pct_Coverage', value = str(coverage)) # only two decimal places.Also float -> str
                etree.SubElement(result,"result_data", type='Depth(min:max:avg)', value =    "{0}:{1}:{2}".format(min_depth, max_depth, avg_depth))
                etree.SubElement(result, "result_data", type='meanQ', value = str(meanQ))
                etree.SubElement(result, "result_data", type='Length', value = str(length))
                etree.SubElement(result, "result_data", type='Coverage_distribution', value = cov_distribution)
                etree.SubElement(result, "result_data", type='Failure', value = failure_tag)
        comment = etree.Comment('(END) Serotype Distinction Results (END)')
        results.append(comment)
        print >> output_file, etree.tostring(root, pretty_print=True)
        output_file.close()

def get_consensus(line):
    """
    Function to provide consensus bp and various quality metrics (depth, freq, insertions, deletions, mixed) for the line provided.

    Arguments:
    line: line from mpileup file ('wzh	13	G	17	.,........,,,,.,,	FFFFFFFFFFFFFFFFB_')
    """
    insertions = re.findall(r'\+\d{1,2}[ACGTNacgtn]+', line.split('\t')[4])
    deletions = re.findall(r'\-\d{1,2}[ACGTacgt]+', line.split('\t')[4]) # this finds pattern -2AA
    mixed = ''
    bps = re.sub(r'\^\S{1}', '', line.split('\t')[4])
    bps = bps.replace('$', '')
    for x in set(insertions):
        no = re.search(r'\d+', x).group() # retrieve number of inserted bps
        nts = re.search(r'[ACGTNacgtn]+', x).group() # retrieve the bps
        new_x = x[:int(no)-len(nts)] if len(nts) != int(no) else x # corrects for cases where +2AAg is selected by reg ex on line 222
        bps = bps.replace(new_x, '')
        if new_x != x:
            for i,y in enumerate(insertions):
                if y==x: insertions[i] = new_x
    for x in set(deletions): bps = bps.replace(x, '') # remove deletion pattern as they do not apply to the current position
    bps = bps.replace('*', '-')
    bps = [x.upper() for x in bps]
    try:
        filtered_bps = [b for i,b in enumerate(bps) if ord(line.split('\t')[5].strip()[i])-33 > 30] # only select bps with mapping quality > 30
    except IndexError:
        pass
    depth = len(filtered_bps)
    if depth == 0:
        consensus = line.split('\t')[2]
        freq = 0
        insertions = deletions = []
    else:
        consensus, cnt = Counter(filtered_bps).most_common()[0] # get the most common bp
        if consensus in [',', '.']:
            consensus = line.split('\t')[2]
            cnt = filtered_bps.count(',')+filtered_bps.count('.')
        freq = round(cnt/float(depth),2)
        if insertions: # adjust consensus frequency to take into account the reads that have an insertion
            # find all insertion with pattern .+2AA or ,+2AA which corresponds to a ref bp with an insertion
            ref_ins = [f for f in re.findall(r'[\.\,]{1}\+\d{1,2}[ACGTNacgtn]+', line.split('\t')[4].upper())] 
            cnt -= len(ref_ins)
            insertions = [f.upper() for f in insertions]
            ins, icnt = Counter(insertions).most_common()[0]
            ins_freq = round(icnt/float(len(line.split('\t')[-1])),2)
            if icnt >= 5 and ins_freq >= 0.75: # if it passes the criteria (d>5, f>0.75) add the insertion to the consensus bp
                consensus += re.search(r'[ACGT]+', ins).group()
                insertions = [ins, ins_freq, icnt]
            elif ins_freq <= 0.75 and ins_freq >= 0.25 and consensus == line.split('\t')[2]: # else recalculate frequence for consensus bp if consensus == ' or .
                freq = round(cnt/float(depth),2)
                mixed = '{nt1}:{cnt1},{nt2}:{cnt2}'.format(nt1=consensus, cnt1=freq, nt2=consensus+re.search(r'[ACGT]+', ins).group(), cnt2=ins_freq)
                if freq < ins_freq:
                    consensus += re.search(r'[ACGT]+', ins).group()
                    insertions = [ins, ins_freq, icnt]
                    freq = ins_freq
        if freq<=0.75 and freq>=0.25 and mixed == '': # report mixed position
            if consensus == line.split('\t')[2] and Counter(filtered_bps).most_common()[1][0] not in [',','.']:
                mixed = '{nt1}:{cnt1},{nt2}:{cnt2}'.format(nt1=consensus, cnt1=freq, nt2=Counter(filtered_bps).most_common()[1][0], cnt2=round(Counter(filtered_bps).most_common()[1][1]/float(depth),2))
                consensus=consensus+'/'+Counter(filtered_bps).most_common()[1][0]
            elif consensus != line.split('\t')[2]:
                second = line.split('\t')[2]
                cnt2 = bps.count(',')+bps.count('.')
                freq2 = round(cnt2/float(len(bps)), 2)
                mixed = '{nt1}:{cnt1},{nt2}:{cnt2}'.format(nt1=consensus, cnt1=freq, nt2=second, cnt2=freq2)
                consensus = consensus+'/'+line.split('\t')[2]

    return consensus.upper(), depth, freq, insertions, deletions, mixed

def get_metrics(all_reads, gene, serotype=None):
    """ Function to calculate various quality metrics for the gene (coverage, avg_depth, min_depth, max_depth, ranges, meanQ, covered_bps, length)

    Arguments:
    all_reads: list with pileup lines for the gene tested
    gene: gene tested
    serotype: serotype tested; optional argument
    """
    if serotype: # only necessary for snps
        length = len(_refseq['{0}_{1}'.format(gene, serotype)]) if '{0}_{1}'.format(gene, serotype) in _refseq.keys() else min([len(_refseq['{0}_{1}'.format(gene,st)]) for st in _serotypes if '{0}_{1}'.format(gene,st) in _refseq.keys()])
    else: # In case of alleles, pseudo or genes the gene name is the key for the _refseq dictionary
        length = len(_refseq[gene])
    coverage = round((len(all_reads)*100)/float(length),1)
    covered_bps = [int(l.split('\t')[1]) for l in all_reads]
    ranges = [] # coverage distribution
    for k, g in groupby(enumerate(covered_bps), lambda (i,x):i-x):
        group_pos = map(itemgetter(1), g)
        ranges.append((group_pos[0], group_pos[-1]))
    depths = [int(l.split('\t')[3]) for l in all_reads]
    avg_depth = round(sum(depths)/float(len(depths)),1) if coverage > 0 else 0
    min_depth = min(depths) if coverage > 0 else 0
    max_depth = max(depths) if coverage > 0 else 0
    avg_qualities = []
    for l in all_reads:
        qualities = [ord(f)-33 for f in l.strip().split('\t')[-1]]
        avg_qualities.append(sum(qualities)/float(len(qualities)))
    meanQ = round(sum(avg_qualities)/float(len(avg_qualities)),1) if avg_qualities != [] else 0

    return coverage, avg_depth, min_depth, max_depth, ranges, meanQ, covered_bps, length

def get_6E_phenotype(snps, reads_all, matches):
    matched={}
    snpsToTest={}
    snpsToTest['wciP']={}
    snpsToTest['wciP'][583]={}
    for st in ['06A', '06B']:
        snpsToTest['wciP'][583][st]=snps['wciP'][583][st]
    snpsToTest['wzy']={}
    snpsToTest['wzy'][657]={'06E': ['ACT', 219, 'T']}
    
    for serotype in ['06A', '06B']:
        matched[serotype]=detect_snps(snpsToTest, reads_all, matches, serotype, 0, 0, [])
    try:
        match = [f for f in matched.keys() if matched[f][0]==1][0]
    except IndexError:
        match= 'n/d'
    snp = [matched[match][2][0]]
    snpsToTest = {}
    snpsToTest['wzy']={}
    snpsToTest['wzy'][657]={'06E': ['ACT', 219, 'T']}
    matched['06E'] = detect_snps(snpsToTest, reads_all, matches, '06E', 0, 0, [])
    if matched['06E'][2] != []: snp.append(matched['06E'][2][0])
    return match, snp

# ============================================================= #
# Pileup functions                                              #
# ============================================================= #
def pileup_investigation(unique, bamfile, reference, serotypes, samtools_path, logger):
    """
    Function to create a mpileup file and test for the various mutations by calling functions detect_genes(), detect_pseudo(),
    detect_alleles() and detect_snps().

    Arguments:
    unique: dictionary containing all mutations the need to be tested to differentiate given serotypes (uploaded from mutationdb.pickle
                     using function parse_serotype_specific_info())
    bamfile: path to sorted bam file created using bowtie_map module
    reference: path to reference fasta file; retrieved using function parse_serotype_specific_info()
    serotypes: list of serotypes to test (genogroup)
    samtools_path: path to samtools
    """
    global _refseq
    global _metrics

    _refseq = parse_reference_fasta(reference)
    _metrics = {} 

    # using samtools mpileup.
    process = subprocess.Popen([samtools_path, 'mpileup', '-B', '-A', '-f', reference, bamfile], stderr=subprocess.PIPE, stdout=subprocess.PIPE)# -A -count anomalous read pairs, -B - disable BAQ computation    and -f FILE - indexed reference sequence file
    reads_all = [l for l in process.stdout] # captures the mpileup output in a list

    # this will select all reads that match the consensus nt at >75% 
    # there is a clause to remove positions with insertions in more than 25% of the reads
    matches = [l for l in reads_all if len(re.findall(r'[\.\,]', l))>len(l.split()[-1])*0.75 and len(re.findall(r'\+{1}\d{1,2}[ACGTacgt]+', l.split('\t')[4]))< len(l.split()[-1])*0.25]     

    sample = os.path.basename(bamfile).split('.')[0]
    detected = {}
    match = []
    matched = {}
    # the unique.keys() indicate what mutation we are testing for ['snps', 'genes', 'alleles' or 'pseudo']
    if 'genes' in unique.keys():
        for serotype in sorted(serotypes):
            if serotype in unique['genes'].keys():
                matched[serotype]=detect_genes(reads_all, serotype, unique['genes'][serotype], serotypes)

    if 'allele' in unique.keys():        
        for serotype in sorted(serotypes):
            mcnt = matched[serotype][0] if serotype in matched.keys() else 0
            tcnt = matched[serotype][1] if serotype in matched.keys() else 0
            unique_matched = matched[serotype][2] if serotype in matched.keys() else []
            matched[serotype] = detect_alleles(unique['allele'], reads_all, serotype, mcnt, tcnt, unique_matched)

    if 'pseudo' in unique.keys():
        for serotype in sorted(serotypes):
            if serotype in unique['pseudo'].keys():
                mcnt = matched[serotype][0] if serotype in matched.keys() else 0
                tcnt = matched[serotype][1] if serotype in matched.keys() else 0
                unique_matched = matched[serotype][2] if serotype in matched.keys() else []
                # check if the gene is a pseudogene
                matched[serotype] = detect_pseudogene(reads_all, matches, unique['pseudo'], serotype, mcnt, tcnt, unique_matched)

    if 'snps' in unique.keys():
        for serotype in sorted(serotypes):
            mcnt = matched[serotype][0] if serotype in matched.keys() else 0
            tcnt = matched[serotype][1] if serotype in matched.keys() else 0
            unique_matched = matched[serotype][2] if serotype in matched.keys() else []
            matched[serotype]= detect_snps(unique['snps'], reads_all, matches, serotype, mcnt, tcnt, unique_matched)

    for serotype in matched.keys(): match.append([serotype, matched[serotype][0], matched[serotype][1], matched[serotype][2]])

    #if 6E add expected phenotype
    scores = [(f[0], f[1]/float(f[2])) for f in match] #creates a list of tuples [(serotype1, score1), (serotype2, score2), ...] where score = matched_hits/total_mutations_tested
    scores = sorted(scores, key=lambda x: x[1], reverse=True)
    if scores[0][0] == '06E':
        phenotype, snp = get_6E_phenotype(unique['snps'], reads_all, matches)
        newHit = '06E({0})'.format(phenotype)
        for i, hit in enumerate(match):
            if hit[0] == '06E':
                replaceWith = [newHit, hit[1]+len([f for f in snp if f!=[]]), hit[2]+len(snp), hit[3]+snp]
                match[i] = replaceWith

    ## export matched into a yaml file to investigate failed samples
    out_yml = open(os.path.join(os.path.dirname(bamfile), 'variant_summary.yml'), 'wb')
    out_yml.write(yaml.dump(match, default_flow_style=True))
    out_yml.close()
    
    return match, sample

def detect_genes(reads_all, serotype, genes, serotypes):
    """
    Function that tests for presence or absence of genes

    Arguments:
    reads_all: list with all pileup lines
    serotype: serotype to test
    genes: dictionary that relates each serotype with expected presence/ absence of a gene {geneA:{serotype1:0, serotype2:1}}
           where 0 = absence and 1 = presence
    serotypes: list of all serotypes tested
    """
    # since this is the first mutation to be tested mcnt, tcnt and unique_matched are defined here
    mcnt = 0
    tcnt =0
    unique_matched = []

    for gene in genes.keys():
        mutation = 'Presence/Absence'
        tcnt += 1
        # get pileup lines for gene
        all_reads = [l for l in reads_all if l.startswith(gene)]
        # get metrics
        coverage, avg_depth, min_depth, max_depth, ranges, meanQ, covered_bps, length = get_metrics(all_reads, gene)
        filtered_coverage = round(len([l for l in all_reads if int(l.split('\t')[3]) >= 5])/float(length), 3)
        _metrics[gene] = [filtered_coverage*100, avg_depth, min_depth, max_depth, ranges, meanQ, length]
        if (filtered_coverage==0 and genes[gene]==0): 
            mcnt+= 1
            unique_matched.append([gene, genes[gene], mutation, 'None', 'genes'])
        else:
            if filtered_coverage >= 0.85 and genes[gene]==1 or filtered_coverage < 0.85 and genes[gene] == 0:
                mcnt+= 1
                unique_matched.append([gene, genes[gene], mutation, 'None', 'genes'])

    return [mcnt, tcnt, unique_matched]

def detect_alleles(genes, reads_all, serotype, mcnt, tcnt, unique_matched):
    """
    Function to determine which allele is present

    Arguments:
    genes: dictionary that relates serotypes to alleles for each gene {geneA:{serotype1:allele1, serotype2:allele2}}
    reads_all: list with all pileup lines
    serotype: serotype to test
    mcnt: number of matched mutations for given serotype
    tcnt: number of all mutations tested for given serotype
    unique_matched: list with all mutations matched for given serotype
    """
    mutation = 'Allele'
    for gene in genes.keys():
        if serotype not in genes[gene].keys(): continue
        failure_tag = 'None'
        # get allele names for this gene
        alleles = [f for f in re.findall(r'[A-Za-z]{3,4}\-{1}\d{1}', ','.join(_refseq.keys())) if f.startswith(gene)]
        tcnt += 1
        coverage = {}
        if gene not in _metrics.keys(): _metrics[gene] = {}
        for allele in alleles:
            # get all pileup line for this allele
            all_reads = [l for l in reads_all if l.startswith(allele)]
            # get metrics
            coverage[allele], avg_depth, min_depth, max_depth, ranges, meanQ, covered_bps, length = get_metrics(all_reads, allele)
            _metrics[gene][allele] = [coverage[allele], avg_depth, min_depth, max_depth, ranges, meanQ, length]
            coverage[allele]=round(coverage[allele]/100, 2)

        if coverage.values().count(1) == 1 and len([f for f in coverage.values() if f >= 0.9])==1: # check for 100% coverage of a single allele
            allele = [f for f in coverage.keys() if coverage[f] == 1][0]
            if genes[gene][serotype] == allele:
                mcnt += 1
                unique_matched.append([gene, allele, coverage, failure_tag, 'allele'])
        elif [f for f in coverage.values() if f >= 0.9] and [f for f in coverage.values() if f < 0.2]: # check for coverage > 90% for one allele and <10% for the second
            allele = [f for f in coverage.keys() if coverage[f] > 0.8][0]
            if genes[gene][serotype] == allele:
                mcnt += 1
                unique_matched.append([gene, allele, coverage, failure_tag, 'allele'])
        elif len([f for f in coverage.values() if f < 0.9]) == len(coverage.values()):
            failure_tag = 'Low coverage: ' + str(coverage)
            unique_matched.append([gene, 'Failed', mutation, failure_tag, 'allele'])
        else: # chimeras or complete mixed sample i.e. 100% coverage for 2 alleles
            failure_tag = 'Mixed: ' + str(coverage)
            unique_matched.append([gene, 'Failed', mutation, failure_tag, 'allele'])

    return [mcnt, tcnt, unique_matched]

def detect_pseudogene(reads_all, matched_reads, pseudogenes, serotype, mcnt, tcnt, unique_matched):
    """
    Function that detects pseudogenes

    Arguments:
    reads_all:list with all pileup lines
    matched_reads: list that contains all pileup lines that match the reference
    pseudogenes: dictionary that relates each serotype with expected presence/ absence of a pseudogene {geneA:{serotype1:0, serotype2:1}}
                 where 0 = functional orf and 1 = pseudogene
    serotype: serotype to test
    mcnt: number of matched mutations for given serotype
    tcnt: number of all mutations tested for given serotype
    unique_matched: list with all mutations matched for given serotype
    """
    for gene in pseudogenes[serotype].keys():
        failure_tag = 'None'
        tcnt += 1

        # get pileup lines for the gene
        all_reads = [l for l in reads_all if l.startswith(gene)]
        matched_all = [l for l in matched_reads if l.startswith(gene)]
        unmatched = set(all_reads)-set(matched_all)

        # ignore positions close -5,+5 bps from a snp; using a sliding window of 10bps with mutation rate>0.2
        ignore_pos = get_ignore_pos(unmatched)

        # get metrics
        coverage, avg_depth, min_depth, max_depth, ranges, meanQ, covered_bps, length = get_metrics(all_reads, gene)
        _metrics[gene] = [coverage, avg_depth, min_depth, max_depth, ranges, meanQ, length]
        coverage = coverage/100
        low_coverage_bps = [int(l.split('\t')[1]) for l in all_reads if int(l.split()[3])<5]

        # positions that are present with depth >= 5
        pass_depth_bps = [l for l in covered_bps if l not in low_coverage_bps]

        # detect pseudogenes
        completed = 'no'
        mutation = ''
        pseudo = ''
        # check if the gene is absent; using 70% rather than 90% as cut off to allow reporting for larger indels if necessary
        if coverage <= 0.7 or (coverage < 0.9 and len(ranges) > 1):
            completed = 'yes'
            if pseudogenes[serotype][gene][0] == 1: # report as cause of pseudogene if serotype tested does not have a functional copy of this gene
                mcnt += 1
                pseudo = 1
                mutation = 'Low coverage'
            else: # or report as failure if functional gene is expected
                failure_tag = 'Low coverage'
                pseudo = 'Failed'
            unique_matched.append([gene, pseudo, 'Pseudogene:'+mutation, failure_tag, 'pseudo'])
        # check if there is 100% match with functional orf
        # used 0.8 cut off instead of 0.9 because of a certain gene that seems to map badly in the beginning and end leading 
        # to wrongful prediction of pseudogene. Confirmed serotype in the lab so I empirically changed the cut off.
        elif len(unmatched) == 0 and coverage >= 0.8: 
            pseudo = 0
            if pseudogenes[serotype][gene][0] == pseudo:
                mcnt += 1
                unique_matched.append([gene, pseudo, 'Pseudogene:'+mutation, failure_tag, 'pseudo'])
            completed = 'yes'
        # look for specific frameshift mutation
        # for example '09A':{'wcjE':[1, [(721,722),'-']]}, '09V':{'wcjE':[0,[(721,722),'G']]}
        elif len(pseudogenes[serotype][gene][1]) == 2:
            completed, matched, mutation, Mixed = check_for_specific_frameshift_mutation(gene, all_reads, matched_all, unmatched, ignore_pos, covered_bps, low_coverage_bps, pseudogenes, serotype)
            if completed == 'yes':
                if matched != '':
                    mcnt += 1
                    unique_matched.append(matched)
            elif completed == 'no' and matched != '':
                if Mixed != []:
                    completed = 'yes'
                    unique_matched.append(matched)
        # look for other ways of inactivating an open reading frame
        if completed == 'no':
            if len(pseudogenes[serotype][gene][1]) == 1:
                # look in specific region
                # i.e. '15B':{'wciZ':[0, [(412,417)]]}, '15C':{'wciZ': [1,[(412,417)]]}
                region = range(pseudogenes[serotype][gene][1][0][0], pseudogenes[serotype][gene][1][0][1])
            else:
                region = []
            if coverage >= 0.9 or (len(ranges) == 1 and coverage > 0.7):
                mutation, pseudo, Mixed = check_for_inactivating_mutations(gene, mutation, region, all_reads, unmatched, low_coverage_bps, ignore_pos, ranges)
                if Mixed != []:
                    mixed_pos = [int(m[0]) for m in Mixed]
                    mixed_pos = sorted(mixed_pos)
                    # identifies larger deletion areas that are mixed. 
                    if mixed_pos[-1]-(mixed_pos[0]-1) == len(mixed_pos):
                        mline = [l for l in all_reads if l.split('\t')[1] == str(mixed_pos[0])][0]
                        consensus, depth, freq, insertions, deletions, mixed = get_consensus(mline)
                        if deletions:
                            if len(set([f.upper() for f in deletions])) == 1:
                                deletion = [f.upper() for f in deletions][0]
                            else:
                                deletion = Counter([f.upper() for f in deletions]).most_common()[0][0]
                            wt = re.search(r'[ACGT]+', deletion).group()
                            dfreq = max([float(re.search(r'-:[0-9\.]{3,4}', m[2]).group().split(':')[1]) for m in Mixed])
                            wt_freq = min([float(re.search(r'[ACGT]{1}:[0-9\.]{3,4}', m[2]).group().split(':')[1]) for m in Mixed])
                            Mixed = [(mixed_pos[0], mixed_pos[-1]), wt+'/'+len(wt)*'-', wt+':'+str(wt_freq), '-'*len(wt)+':'+str(dfreq)]
                    failure_tag = 'Mixed:'+ ','.join([str(m) for m in Mixed])
                    pseudo = 'Failed'
                    if mutation == []: mutation = 'None'
                    unique_matched.append([gene, pseudo, 'Pseudogene:'+str(mutation), failure_tag, 'pseudo'])
            if pseudogenes[serotype][gene][0] == pseudo:
                mcnt += 1
                if mutation == []: mutation = 'None'
                unique_matched.append([gene, pseudo, 'Pseudogene:'+str(mutation), failure_tag, 'pseudo'])

    return [mcnt, tcnt, unique_matched]

def detect_snps(snps, reads_all, matches, serotype, mcnt, tcnt, unique_matched):
    """
    Function that detect for snps in the pileup file.

    Arguments:
    snps: dictionary that contains all snps; format: {geneA: {pos1: {serotype1:[codon, aa_pos, aa], serotype2:...}, pos2: {...}}, geneB: {..}...}
    reads_all: list with all pileup lines
    matches: list that contains all pileup lines that match the reference 
    serotype: serotype tested
    mcnt: number of matched mutations for given serotype
    tcnt: number of all mutations tested for given serotype
    unique_matched: list with all mutations matched for given serotype
    """
    for gene in snps.keys():
        failure_tag = 'None'
        mutation = 'snps'
        serotypes = []
        
        if gene not in _metrics.keys(): _metrics[gene] = {}
        if serotype not in [serotypes+a.keys() for a in snps[gene].values()][0]:
            continue
        # select only pileup lines for the gene in question
        all_reads = [l for l in reads_all if l.startswith('{0}_{1}'.format(gene, serotype))]
        matched_reads = [l for l in matches if l.startswith('{0}_{1}'.format(gene, serotype))]
        # get metrics
        coverage, avg_depth, min_depth, max_depth, ranges, meanQ, covered_bps, length = get_metrics(all_reads, gene, serotype)
        _metrics[gene][serotype] = [coverage, avg_depth, min_depth, max_depth, ranges, meanQ, length] # for snps each serotype has a distinct gene sequence
        coverage=round(coverage/100, 2)
        # find all positions with depth < 5
        low_coverage_bps = [int(l.split('\t')[1]) for l in all_reads if int(l.split()[3])<5]
        if coverage < 0.9 and serotype not in ['11A', '11D', '11F']: # coverage threshold for snps
            failure_tag = 'Low Coverage'
            unique_matched.append([gene, 'Failed', mutation, failure_tag, 'snps'])
            for pos in snps[gene].keys():
                if serotype in snps[gene][pos].keys(): tcnt += 1
            continue
        for pos in snps[gene].keys():
            failure_tag = 'None'
            if serotype not in snps[gene][pos].keys():
                continue
            elif pos not in covered_bps:
                failure_tag = 'Low Coverage'
                unique_matched.append([gene, 'Failed', mutation, failure_tag, 'snps'])
                tcnt += 1
                continue
            tcnt += 1
            aapos = int(pos)/3
            # checks that all nts within the codon for this aa are in positions with depth > 5
            if len([f for f in range(int(aapos)*3+1,int(aapos)*3+1+3) if f in covered_bps and f not in low_coverage_bps]) == 3:
                if len(all_reads) == len(matched_reads): # if all reasds match the reference then the gene sequence corresponds to the serotype tested
                    mcnt += 1
                    unique_matched.append([gene, pos, snps[gene][pos][serotype], failure_tag, 'snps'])
                else: # search for specific position
                    if type(pos) == int: # snps
                        aapos = pos/3
                        expectedCodon = snps[gene][pos][serotype][0]
                        #get pileup lines for specific codon -> pos in range(int(aapos)*3+1,int(aapos)*3+1+3)
                        rec_all = [l for l in all_reads if int(l.split('\t')[1]) in range(int(aapos)*3+1,int(aapos)*3+1+3)]
                        rec_matched = [l for l in matched_reads if int(l.split('\t')[1]) in range(int(aapos)*3+1,int(aapos)*3+1+3)]
                        rec_mixed = [l for l in rec_all if len(re.findall(r'[\.\,]', l.split('\t')[4]))<len(l.split()[-1])*0.75 and len(re.findall(r'[\.\,]', l.split('\t')[4]))>len(l.split()[-1])*0.25]
                        seq = ''
                        for l in rec_all:
                            consensus, depth, freq, insertions, deletions, mixed = get_consensus(l)
                            if freq > 0.75 and consensus in ['A', 'C', 'G', 'T', '-'] and depth>=5:
                                seq += consensus
                        if len(rec_all) == len(rec_matched) == 3 and rec_mixed == []: # codon matched the reference
                            mcnt += 1
                            unique_matched.append([gene, pos, snps[gene][pos][serotype], failure_tag, 'snps'])
                        elif expectedCodon == seq:
                            mcnt += 1
                            unique_matched.append([gene, pos, snps[gene][pos][serotype], failure_tag, 'snps'])
                        else:
                            # check if the position is mixed
                            if rec_mixed:
                                mx_sum = []
                                for l in rec_mixed:
                                    GT = {}
                                    for nt in ['A','C','G','T']:
                                        if nt == l.split('\t')[2]:
                                            GT[nt] = len(re.findall(r'[\.\,]', l.split('\t')[4]))
                                        else:
                                            GT[nt] = len(re.findall(r'[{0}{1}]'.format(nt, nt.lower()), l.split('\t')[4]))
                                    if [f for f in GT.values() if f!=max(GT.values()) and f/float(max(GT.values()))>0.25]:
                                        mx_sum.append(str((l.split('\t')[1], GT)))
                                failure_tag = 'Mixed: ' + ';'.join(mx_sum) if mx_sum!= [] else ''
                                unique_matched.append([gene, pos, 'Failed', failure_tag, 'snps'])
            else: # report failure for specific position
                if len([f for f in range(int(aapos)*3+1,int(aapos)*3+1+3) if f in covered_bps]) != 3:
                    failure_tag = 'SNP position {0} not covered'.format(pos)
                elif len([f for f in range(int(aapos)*3+1,int(aapos)*3+1+3) if f not in low_coverage_bps]) != 3:
                    failure_tag = 'SNP position {0}: Depth below 5'.format(pos)
                unique_matched.append([gene, pos, 'Failed', failure_tag, 'snps'])

    return [mcnt, tcnt, unique_matched]

# ============================================================= #
# Small functions for detect_pseudogene()                       #
# ============================================================= #

def check_for_inactivating_mutations(gene, mutation, region, all_reads, unmatched, low_coverage_bps, ignore_pos, ranges):
    """
    Function that determines whether a gene is functional or not, looking for a range of
    inactivating mechanisms (frameshift mutations, early stop codons, truncations) throughout
    the whole length of the gene sequence.

    Arguments:
    gene: gene tested
    mutation: a list of mutations (mixed position tested previously with check_for_specific_frameshift_mutation())
              or a string corresponding to a single mutation that failed for other reasons (coverage, mutation rate)
    region: either an empty list or range (list) of positions if provided by the database
    all_reads: list of pileup lines for the gene tested
    unmatched: list of pileup lines that do not match the reference
    ignore_pos: positions with more than 2 snps within a window of 10bps; mutation rate > 0.2
    low_coverage_bps: positions with depth < 5
    ranges: coverage distribution for gene tested
    """
    frameshift = ''
    del_pos = []
    seq = ''
    Mixed = []
    if type(mutation) != list: mutation = [] 
    start = int(all_reads[0].split('\t')[1]) - 1
    end = int(all_reads[-1].split('\t')[1]) - 1
    indels = get_indels(unmatched)
    # deals with cases where there is problem mapping at the start and end of the gene sequence.
    if start >= 1 and start % 3 != 0:
        for i in range(1,3):
            # removes 1-2 bases to make sure that when extracting the seq it would be in frame.
            if (int(start)+i)%3 == 0:
                all_reads = all_reads[i:]
    for i,l in enumerate(all_reads):
        if int(l.split('\t')[1]) >= end:
            continue
        if int(l.split('\t')[1]) in low_coverage_bps:
            seq += l.split('\t')[2]
            continue
        consensus, depth, freq, insertions, deletions, mixed = get_consensus(l)
        if int(l.split('\t')[1]) in ignore_pos:
            # if the snp is at the beginning or end of the covered area
            # ingore it as the ends of the reads seem to be variable as a byproduct
            # of sequencing. Discuss with Anthony??? Example COPENHAGEN-33A-2 and wcjE
            seq+=l.split('\t')[2]
        elif freq >= 0.75 and consensus in ['A', 'C', 'G', 'T'] and depth>= 5:
            seq+= consensus.upper()
            if l in unmatched and consensus != l.split('\t')[2]:
                mutation.append('{0}{1}{2}'.format(l.split('\t')[2],l.split('\t')[1],consensus))
        elif freq >= 0.75 and consensus == '-' and depth>= 5:
            mutation.append(['{0}del{1}'.format(l.split('\t')[1], l.split('\t')[2])]) 
            continue
        elif len(consensus) > 1 and freq >= 0.75 and mixed == '':
            seq+= consensus.upper()
            inserted_bps = re.sub('\+\d{1,2}', '', insertions[0])
            mutation.append('{0}ins{1}'.format(l.split('\t')[1], inserted_bps))
            if len(inserted_bps)%3 != 0:
                frameshift = '{0}ins{1}'.format(l.split('\t')[1], inserted_bps)
        else:
            seq+= l.split('\t')[2]
        if mixed and int(l.split('\t')[1]) in region and int(l.split('\t')[3]) > 5:
            Mixed.append([l.split('\t')[1],consensus, mixed])
    # translate the consensus sequence. Adjust according to ranges to remain on the same orf
    ### removed since at lines 746-750 the script removes starting bases out of frame to make sure the seq is within the correct ORF
    #if ranges[0][0]%3 == 0: protein = Seq(seq[ranges[0][0]%3+1:],generic_dna).translate()
    #elif ranges[0][0]%3 == 1: protein = Seq(seq[ranges[0][0]%3-1:],generic_dna).translate()
    #elif ranges[0][0]%3 == 2: protein = Seq(seq[ranges[0][0]%3:],generic_dna).translate()
    protein = Seq(seq,generic_dna).translate()
    orfs = protein.split('*')
    orfs.sort(key=lambda x:len(x), reverse=True)
    # define the expected dna and protein sequene for the functional gene
    functional_orf = _refseq[gene].translate()
    functional_seq = _refseq[gene]
    serotype_psd = 0
    # compare the consensus sequence to the expected functional sequence
    if len(orfs[0]) in range(len(functional_orf)-15,len(functional_orf)+1) and (str(orfs[0]) in str(functional_orf) or str(orfs[0][orfs[0].find('M'):]) in str(functional_orf)):
        pseudo = 0
    elif str(functional_orf).split('*')[0] in protein:
        pseudo = 0 
    elif seq in functional_seq and Mixed == []:
        pseudo = serotype_psd
    elif protein.count('*') <= 1:
        alignments = pairwise2.align.globalxx(str(orfs[0]), str(functional_orf).split('*')[0])
        identity = alignments[0][2]/float(alignments[0][4])
        pseudo = 0 if identity >= 0.9 else 1
    else:
        pseudo = 1
        if frameshift == '' and indels!=[]: # if indels are present and no frameshift has been reported yet then it must be deletions
            first_stop_cdn_pos = protein.find('*')
            del_lines = [f for f in indels if int(f.split('\t')[1]) < first_stop_cdn_pos*3]
            del_pos = [int(f.split('\t')[1]) for f in del_lines if f.split('\t')[4].count('*') > len(f.split('\t')[5].strip())*0.75]
            del_pos = sorted(del_pos)
            if del_pos and del_pos[-1]-(del_pos[0]-1) == len(del_pos):
                del_line = [l for l in all_reads if l.split('\t')[1] == str(del_pos[0]-1)][0]
                if re.findall(r'[\-]{1}\d{1,2}[acgtACGT]+', del_line.split('\t')[4]) == []:
                    # deals with cases where part of the deletion falls sort of the 0.75 cut off so reducing it to 0.70 to check
                    del_pos = [int(f.split('\t')[1]) for f in del_lines if f.split('\t')[4].count('*') > len(f.split('\t')[5].strip())*0.70]
                    del_pos = sorted(del_pos)
                    if del_pos[-1]-(del_pos[0]-1) == len(del_pos): del_line = [l for l in all_reads if l.split('\t')[1] == str(del_pos[0]-1)][0]
                consensus, depth, freq, insertions, deletions, mixed = get_consensus(del_line)
                try:
                    deletion = [f.upper() for f in deletions][0] if len(set([f.upper() for f in deletions])) == 1 else Counter([f.upper() for f in deletions]).most_common()[0][0]
                    wt = re.search(r'[ACGT]+', deletion).group()
                    frameshift = '{0}del{1}'.format(del_pos[0], wt)
                except IndexError:
                    pass
            if frameshift != '': mutation = frameshift
        elif len(ranges)> 1:
            mutation = 'Truncation'
        elif frameshift != '':
            mutation = frameshift
        elif frameshift == '' and indels==[] and mutation != []:
            snp_pos = [int(f[1:-1]) for f in mutation]
            stop_pos = protein.find('*')
            if str(Seq(seq[stop_pos*3:stop_pos*3+3], generic_dna).translate()) == '*':
                wt_aa = str(functional_seq[stop_pos*3:stop_pos*3+3].translate())
                mutation = '{0}{1}*'.format(wt_aa, stop_pos)

    return mutation, pseudo, Mixed

def check_for_specific_frameshift_mutation(gene, all_reads, matched_all, unmatched, ignore_pos, covered_bps, low_coverage_bps, pseudogenes, serotype):
    """
    Function that checks for a specific frameshift mutation known to be responsible for inactivating a
    functional orf in serotype tested
    
    Arguments:
    all_reads: list of pileup lines for the gene tested
    matched_all: list of pileup lines that match the reference
    unmatched: list of pileup lines that do not match the reference
    ignore_pos: positions with more than 2 snps within a window of 10bps; mutation rate > 0.2
    covered_bps: positions covered in the mpileup file
    low_coverage_bps: positions with depth < 5
    pseudogenes: dictionary that relates each serotype with expected presence/ absence of a pseudogene 
                 i.e. '09A':{'wcjE':[1, [(721,722),'-']]}, '09V':{'wcjE':[0,[(721,722),'G']]}
                 where 0 = functional orf and 1 = pseudogene
    serotype: serotype tested
    """
    # get info for serotype tested
    region, mut = pseudogenes[serotype][gene][1]
    start,end = region
    pseudo_status = pseudogenes[serotype][gene][0]
    completed = 'yes'
    Mixed = []
    failure_tag = ''
    
    # define frameshift mutation provided by the database
    frameshift = [pseudogenes[st][gene][1] for st in pseudogenes.keys() if pseudogenes[st][gene][0] == 1][0][1]
    wt = [pseudogenes[st][gene][1] for st in pseudogenes.keys() if pseudogenes[st][gene][0] == 0][0][1]
    if frameshift.count('-') == 0 and len(wt) != len(frameshift):
        mutation = '{0}ins{1}'.format(start+1, frameshift[1:])
    elif frameshift.count('-') != 0:
        mutation = '{0}del{1}'.format(start+1, wt)
    else:
        mutation = '{0}{1}*'.format(Seq(wt, generic_dna).translate(), start/3)

    # check for failures in given region
    if [f for f in range(start+1, end+1) if f in ignore_pos]:
        completed = 'no'
        failure_tag = 'Frameshift mutation {0} in high variability site'.format(mutation)
        # Sliding window 10bps - mutation rate>0.2
    elif [f for f in range(start+1, end+1) if f not in covered_bps]: #not covered
        completed = 'no' 
        failure_tag = 'Frameshift mutation {0} not covered'.format(mutation)
    elif [f for f in range(start+1, end+1) if f in low_coverage_bps]: # depth<5
        completed = 'no'
        failure_tag = 'Frameshift mutation {0} in low depth area'.format(mutation)

    # if the region specified passed all criteria continue with detection
    matched = ''
    if completed != 'no': 
        rec_all = [l for l in all_reads if int(l.split('\t')[1]) in range(int(start)+1,int(end)+1)]
        rec_matched = [l for l in matched_all if int(l.split('\t')[1]) in range(int(start)+1,int(end)+1)]
        if len(rec_all) == len(rec_matched) == len(range(start, end)):
            pseudo = 0
            if pseudogenes[serotype][gene][0] == pseudo and unmatched == []:
                matched = [gene, pseudo, 'Pseudogene: None', failure_tag, 'pseudo']
            else:
                completed = 'no'
        elif len(rec_all) != len(rec_matched):
            if re.findall(r'[ACGT]+', mut) == [] and mut.count('-') == len([f for f in rec_all if f.split('\t')[4].count('*')>len(f.split()[-1])*0.75]) and pseudogenes[serotype][gene][0] == 1:
                # frameshift detected 
                pseudo = 1
                matched = [gene, pseudo, 'Pseudogene:'+mutation, failure_tag, 'pseudo']
            elif re.findall(r'[ACGT]+', mut) == [] and mut.count('-') == len([f for f in rec_all if f.split('\t')[4].count('*')>=len(f.split()[-1])*0.2]):
                # frameshift detected but in mixed state
                seq = ''
                for l in rec_all:
                    consensus, depth, freq, insertions, deletions, mixed = get_consensus(l)
                    if mixed: Mixed.append([l.split('\t')[1], consensus, mixed])
                    if freq > 0.75 and consensus in ['A', 'C', 'G', 'T', '-'] and depth>=5:
                        seq+= consensus
                    elif freq<=0.75 and freq>=0.25 and depth>=5:
                        seq+=consensus
                    elif depth < 5:
                        completed ='no'
                        Mixed.append([l.split('\t')[1],'Low coverage: Depth below 5'])
                mutation = [[mutation, Mixed]]
            elif re.findall(r'[ACGT]+', mut) != []:
                mutations = [pseudogenes[st][gene][1][1] for st in pseudogenes.keys() if pseudogenes[st][gene][1] != [] ]
                seq = ''
                for l in rec_all:
                    consensus, depth, freq, insertions, deletions, mixed = get_consensus(l)
                    if depth < 5:
                        pseudo = 'Failed'
                        completed = 'no'
                        matched = [gene, pseudo, 'Pseudogene', mutation, 'Depth below 5', 'pseudo']
                        break
                    if freq > 0.75 and consensus in ['A', 'C', 'G', 'T', '-']:
                        seq+= consensus
                    elif freq<=0.75 and freq>=0.25:
                        if mixed: Mixed.append([l.split('\t')[1],consensus, mixed])
                        seq+=consensus
                        completed = 'no'
                    #elif freq > 0.75 and consensus in mutations:
                    #    seq+=consensus
                if seq == mut and Mixed == []:
                    pseudo = pseudogenes[serotype][gene][0]
                    matched = [gene, pseudo, 'Pseudogene:'+mutation, failure_tag, 'pseudo']
                elif seq in mutations and Mixed != []:
                    pseudo = 'Failed'
                    failure_tag = 'Mixed:'+ ','.join([str(m) for m in Mixed])
                    matched = [gene, pseudo, 'Pseudogene:'+mutation, failure_tag, 'pseudo']
                    completed = 'no'
                elif seq not in mutations:
                    completed ='no'

    return completed, matched, mutation, Mixed

def get_indels(unmatched_lines):
    """
    Small function that returns the pileup lines with indels seen in more than 75% of the reads at each position

    Arguments:
    unmatched_reads: list of pileup lines for positions that do not match the reference
    """
    indels = []
    for f in unmatched_lines:
        ins = re.findall(r'[\+]{1}\d{1,2}[acgtACGT]+', f.split('\t')[4])
        expected_threshold = len(f.split('\t')[5].strip())*0.75 # 75% of the reads
        if int(f.split('\t')[3]) > 5 and (len(ins) > expected_threshold or f.split('\t')[4].count('*') > expected_threshold):
            indels.append(f)

    return indels

def get_ignore_pos(unmatched_lines):
    """
    Small function that returns snp positions within close proximity to other snps (sliding window of 10 bps)

    Arguments:
    unmatched_lines: list of pileup lines for positions that do not match the reference
    """
    ignore_pos = []
    unmatched_pos = [int(f.split('\t')[1]) for f in unmatched_lines]
    for i,f in enumerate(unmatched_pos):
        sliding_window = [m for m in range(f-5,f+5) if m in unmatched_pos]
        mutation_rate = len(sliding_window)/len(range(f-5,f+5))
        line = [g for g in unmatched_lines][i]
        if mutation_rate >0.2 and line.split('\t')[4].count('*') == 0:
            ignore_pos.append(f)

    return ignore_pos
  

# ============================================================= #
# Main function                                                #
# ============================================================= #

def find_serotype(fastq_path, serotypes, reference_directory, output_dir, bowtie_path, samtools_path, clean_bam, logger, workflow='', version = ''):
    """
    This is the main function that calls other functions to determine serotype.

    Arguments:
    fastq_path: path fastq file; either input directory containing two fastq files or the filepath to the 1.fastq.gz (or 1.fastq) 
    serotypes: list of closely related serotypes in this format [07A, 07F, 40] i.e. all serotypes with numbers and letters need to follow this format \d{2}[A-Z]{1}
    output_dir: path to output directory 
    bowtie: path to bowtie
    samtools: path to samtools
    """
    # make output_directory
    outdir = os.path.join(output_dir, 'SNP_based_serotyping')
    if not os.path.exists(outdir): os.makedirs(outdir)

    # extract path to fastq file; two possibilities either already provided or input directory provided
    glob_pattern = "*fastq*"
    if os.path.isdir(fastq_path): # if input directory provided
        input_dir = fastq_path
        fastq_files = glob.glob(os.path.join(fastq_path, glob_pattern))
        if len(fastq_files) != 2 and fastq_files != []:
            print "Unexpected number (" + str(len(fastq_files)) + ") of fastq files"
            sys.exit(1)
        if fastq_files == []:
            fastq_files = glob.glob(os.path.join(fastq_path, '*.fastq*'))
            if fastq_files == []:
                print "No fastq files in the directory provided:", fastq_path
                sys.exit(1)
        fastq = [f for f in fastq_files if f.find('1.fastq') != -1][0]
    else:
        input_dir = os.path.dirname(fastq_path)
        if input_dir == '': input_dir = os.getcwd()
        if fastq_path.endswith('fastq.gz') or fastq_path.endswith('fastq'):
            fastq = fastq_path
        else:
            print 'Please provide a fastq file'
            sys.exit(1)

    # read the reference files (reference fasta and mutationdb.pickle)associated with given serotypes
    unique, reference_fasta_file, serotypes = try_and_except(output_dir+"/logs/SNP_based_serotyping.stderr", parse_serotype_specific_info, serotypes, reference_directory, outdir, fastq, logger, workflow, version)
    if unique != None:
        sorted_bamfile, reference_fasta_file = try_and_except(input_dir+"/logs/SNP_based_serotyping.stderr", bowtie_map.mapping, fastq, reference_fasta_file, bowtie_path, samtools_path, outdir, logger)
        match, sample = try_and_except(output_dir+"/logs/SNP_based_serotyping.stderr", pileup_investigation, unique, sorted_bamfile, reference_fasta_file, serotypes, samtools_path, logger)
        try_and_except(output_dir+"/logs/SNP_based_serotyping.stderr", output_xml, match, outdir, sample, unique, logger, workflow, version)
        if clean_bam == True: # remove bam files if -c option is selected
            bamfiles = glob.glob(sorted_bamfile+'*')
            for file in bamfiles: os.remove(file)
    # clean up all unnecessary files
    if os.path.exists(outdir):
        tmp = glob.glob(os.path.join(outdir, '*tmp*'))
        if tmp:
            clean_up(tmp[0], logger)

  
