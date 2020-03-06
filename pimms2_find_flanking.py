# from argparse import _SubParsersAction
# import configparser
# import pathlib
from fuzzysearch import find_near_matches
import os
import time
import sys
import re
import configargparse
import datetime
import multiprocessing
import glob
import gzip
import pysam
import subprocess

import logging
import logging.handlers


# from time import sleep
# from random import random, randint

# from pathlib import Path


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def deleteFileList(fileList):
    for filePath in fileList:
        try:
            os.remove(filePath)
        except OSError:
            print("Error while deleting file")


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise configargparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def concat_fastq(flanking_fastq_list, label, fq_file_suffix):
    concat_fastq_result_filename = label + '_RX_concat' + fq_file_suffix
    print(concat_fastq_result_filename)
    print(" ".join(flanking_fastq_list))

    with pysam.FastxFile(flanking_fastq_list[0]) as fin, gzip.open(concat_fastq_result_filename, mode='wt') as fout:
        for entry in fin:
            fout.write((str(entry) + '\n'))

    flanking_fastq_list.pop(0)
    for result_fq in flanking_fastq_list:
        print('FQRESULT' + result_fq + '\n')
        with pysam.FastxFile(result_fq) as fin, gzip.open(concat_fastq_result_filename, mode='at') as fout:
            for entry in fin:
                fout.write((str(entry) + '\n'))

    return (concat_fastq_result_filename)


def run_minimap2(flanking_fastq_concat_result, sam_output_result, genome_fasta):
    stream = os.popen('minimap2 --version')
    output = stream.read()
    print(output)
    # process = subprocess.Popen(['minimap2', '--version'],
    print(' '.join(['minimap2', '-x', 'sr', '-a', '-o', sam_output_result, genome_fasta, flanking_fastq_concat_result,
                    '--secondary=no', '--sam-hit-only']))
    process = subprocess.Popen(
        ['minimap2', '-x', 'sr', '-a', '-o', sam_output_result, flanking_fastq_concat_result, '--secondary=no',
         '--sam-hit-only'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stdout.decode('utf-8'))
    print(stderr.decode('utf-8'))


pimms_mls = """
===========================================================================================================
Pragmatic Insertional Mutation Mapping system (PIMMS) mapping pipeline 2
===========================================================================================================
       o         o
       o         o
      //        //
     //        //
   |_||_|    |_||_|   @@@@@  @@@@@@  @@     @@  @@     @@   @@@@@@    @@@@@@
   |@||@|    |@||@|   @@  @@   @@    @@@@ @@@@  @@@@ @@@@  @@        @@    @@
   |@||@|    |@||@|   @@@@@    @@    @@ @@@ @@  @@ @@@ @@   @@@           @@
   |@||@|    |@||@|   @@       @@    @@  @  @@  @@  @  @@     @@@        @@
   |@@@@|    |@@@@|   @@       @@    @@     @@  @@     @@       @@     @@
   |@@@@|    |@@@@|   @@     @@@@@@  @@     @@  @@     @@  @@@@@@@   @@@@@@@@
===========================================================================================================
python PIMMS2 test script.....
===========================================================================================================\n"""

# Construct the argument parser
# ap = argparse.ArgumentParser(description='PIMMS2 fastq sequence processing', prog="demo_pimms2",
#                              epilog='\n\n*** N.B. This is a development version ***\n \n')
# modes: _SubParsersAction = ap.add_subparsers()

ap = configargparse.ArgumentParser(  # description='PIMMS2 sam/bam processing',
    prog="demo_pimms2_find_flank",
    add_config_file_help=False,
    config_file_parser_class=configargparse.DefaultConfigFileParser,
    epilog="\n\n*** N.B. This is a development version ***\n \n ",
    description='''description here'''
)
ap.add_argument('-v', '--version', action='version', version='%(prog)s 2.0.1 demo')
modes = ap.add_subparsers(parser_class=configargparse.ArgParser)

# modes.required = False
findflank = modes.add_parser("find_flank", add_config_file_help=False,
                             help="Mode: extract insertion site coordinates from sam file",
                             description="Args that start with '--' (eg. --sam) can also be set in a config file (specified via -c)")
samcoords = modes.add_parser("sam_extract", help="Mode: extract insertion site coordinates from sam file")
otherstuff = modes.add_parser("other_stuff", help='Mode: do other good PIMMS related stuff')
# Add the arguments to the parser

findflank.add_argument("-c", "--config", required=False, is_config_file=True,  # dest='config_file',
                       # type=str, default='',
                       metavar='pimms2.config',
                       help="use parameters from config file")
# findflank.add_argument("--single", required=False, action='store_true',
#                       help="illumina single end reads (un-paired) [default: illumina paired end]")
# findflank.add_argument("--rmvector", required=False, action='store_true',
#                       help="attempt removal of contaminating vector")
findflank.add_argument("--nano", required=False, action='store_true', default=False,
                       help="override with settings more suitable for nanopore")
findflank.add_argument("--fasta", required=False, nargs=1, metavar='ref_genome.fasta', type=extant_file,
                       help="fast file for reference genome ")
findflank.add_argument("--nomap", required=False, action='store_true', default=False,
                       help="do not run mapping step")
findflank.add_argument("--rmfiles", required=False, action='store_true', default=False,
                       help="remove intermediate files")
findflank.add_argument("--lev", required=False, action='store_true', default=False,
                       help="use Levenshtein distance of 1")
findflank.add_argument("--sub", required=False, nargs=1, type=int, default=1,
                       help="number of permitted base substitutions in motif match [1]")
findflank.add_argument("--insert", required=False, nargs=1, type=int, default=0,
                       help="number of permitted base insertions in motif match [0]")
findflank.add_argument("--del", required=False, nargs=1, type=int, default=0, dest='deletion',
                       help="number of permitted base insertions in motif match [0]")
findflank.add_argument("--in_dir", required=True, nargs=1, dest='in_dir', type=extant_file,
                       help="directory containing input fastq files (assumed to match '*q.gz' or '*.fastq')")
findflank.add_argument("--fwdrev", required=False, nargs=1, type=str, default=['_R1_,_R2_'],
                       help="text substring to uniquely identify illumina fwd/rev paired fastq files ['_R1_,_R2_']")
findflank.add_argument("--out_dir", required=False, nargs=1, metavar='DIR', default=[''],  # dest='out_dir'
                       action='store',
                       # (['pimms2_' + time.strftime("%y%m%d_%H%M%S")]),
                       help="directory to contain fastq files results")
findflank.add_argument("--prefix", required=False, nargs=1, type=str, default="pimms2_condition",
                       help="prefix for output files")
findflank.add_argument("--cpus", required=False, nargs=1, type=int, default=int(os.cpu_count() / 2),
                       help="number of processors to use [(os.cpu_count() / 2)] ")
findflank.add_argument("--max", required=True, nargs=1, type=int, default=60,
                       help="clip results to this length [illumina:60/nano:100")
findflank.add_argument("--min", required=True, nargs=1, type=int, default=25,
                       help="minimum read length [illumina:60/nano:100")
findflank.add_argument("--motif1", required=False, nargs=1, type=str, default=['TCAGAAAACTTTGCAACAGAACC'],
                       help="IS end reference motif1 [pGh9:TCAGAAAACTTTGCAACAGAACC]")
findflank.add_argument("--motif2", required=False, nargs=1, type=str, default=['GGTTCTGTTGCAAAGTTTAAAAA'],
                       help="IS end reference motif2 [pGh9:GGTTCTGTTGCAAAGTTTAAAAA]")

findflank.add_argument("--label", required=False, nargs=1, metavar='condition_name', default=[''],
                       help="text tag to add to results file")

parsed_args = ap.parse_known_args()

# exit and print sort help message if no mode/arguments supplied
if len(sys.argv) <= 2:
    ap.print_usage()
    sys.exit(1)

# do command line processing


print(parsed_args)
print("----------")
# print(ap.format_help())
print("----------")
print(ap.format_values())  # useful for logging where different settings came from

# print((vars(parsed_args)))

# config_file = parsed_args.config_file[0]

# construct config parser
# p2config = configparser.ConfigParser()

if not parsed_args[0].nomap:
    if parsed_args[0].fasta is None:
        ap.error("unless the --nomap flag is used please supply a sequence file e.g:  --fasta contigs.fasta")
    elif not parsed_args[0].label:
        ap.error("unless the --nomap flag is used please supply a text label string  --label cond_01")
    else:
        print("refseq provided: " + parsed_args[0].fasta[0])


else:
    print('Skipping sequence mapping ')

if parsed_args[0].out_dir[0]:
    out_dir = parsed_args[0].out_dir[0]
# print('\ncreating result dir: ' + out_dir + '\n')
else:
    out_dir = 'pimms2_' + time.strftime("%y%m%d_%H%M%S")
# print('\ncreating result dir: ' + out_dir + '\n')
# createFolder(out_dir)

if os.path.isdir(out_dir):
    print('\nresult dir exists\n')
else:
    print('\ncreating result dir: ' + out_dir + '\n')
createFolder(out_dir)

fwdrev_wc = parsed_args[0].fwdrev[0].strip("'\"").split(',')

# exit(0)

# print(pimms_mls)

ncpus = int(parsed_args[0].cpus)
nano = parsed_args[0].nano
# experimental decontaminate transposon/vector sequence
# not currently effective try another implementation when time allows?
decontam_tranposon = False
fuzzy_levenshtein = parsed_args[0].lev

# if decontam_tranposon == False:
#     decon_tag = "nodecon"
# else:
#     decon_tag = "decon"

# fq_result_suffix = "_pimmsout_trim100_nodecon.fastq"

if nano:  # nano == True
    fuzzy_levenshtein = True
    l_dist = 2  # maximum Levenshtein Distance
    min_length = 30
    max_length = 150
    print('overriding with Nanopore appropriate settings: Levenshtein distance of ' + str(
        l_dist) + ' + sequence length min = ' + str(min_length) + ', max = ' + str(max_length))
else:
    subs = parsed_args[0].sub[0]
    l_dist = 1  # maximum Levenshtein Distance
    insrt = parsed_args[0].insert[0]
    dels = parsed_args[0].deletion[0]
    min_length = parsed_args[0].min[0]
    max_length = parsed_args[0].max[0]

if fuzzy_levenshtein:
    fq_result_suffix = ("_pimms2out_trim" + str(max_length) + "_lev" + str(l_dist) + ".fastq.gz")
elif insrt > 0 | dels > 0:
    fq_result_suffix = ("_pimms2out_trim" + str(max_length) + "_sub" + str(subs) + "_ins" + str(insrt) + "_del" + str(
        dels) + ".fastq.gz")
else:
    fq_result_suffix = ("_pimms2out_trim" + str(max_length) + "_sub" + str(subs) + ".fastq.gz")

sam_result_suffix = re.sub('.fastq.gz', '.sam', fq_result_suffix)

trans = str.maketrans('ATGCN', 'TACGN')  # complement DNA lookup

# contam1 = 'GATGCTCTAGAGCATTCTCT'
# contam2 = 'CATTCTCTCCATCAAGCTAT'
# contam1rc = contam1.translate(trans)[::-1]  # reverse complement ([::-1] -> reverse)
# contam2rc = contam2.translate(trans)[::-1]

# qry1 = "TCAGAAAACTTTGCAACAGAACC"
# qry2 = "GGTTCTGTTGCAAAGTTTAAAAA"
qry1 = parsed_args[0].motif1[0].strip("'\"")
qry2 = parsed_args[0].motif2[0].strip("'\"")
print(str(qry1))
print(str(qry2))
# revcomp using maketrans lookup and a string reverse
qry1rc = qry1.translate(trans)[::-1]  # reverse complement transposon motif1 ([::-1] -> reverse)
qry2rc = qry2.translate(trans)[::-1]  # reverse complement transposon motif2

max_length_index = max_length - 1

fastq_dir = os.path.join(parsed_args[0].in_dir[0], '')


#

def pimms_fastq(fq_filename, fqout_filename):
    print(fq_filename + "\n" + fqout_filename + "\n")
    count = 0
    countq1 = 0
    countq2 = 0
    countq1q2 = 0
    countq1rc = 0
    countq2rc = 0
    countq1rcq2rc = 0
    hit_but_short_q1_q2 = 0
    hit_q1_q2 = 0
    countq2rcq1rc = 0
    hit_but_short_q2rc_q1rc = 0
    hit_q2rc_q1rc = 0
    wrongq2q1 = 0
    wrongq1rcq2rc = 0
    countqqrc = 0
    countqmulti = 0
    hit_but_short_q1_only = 0
    hit_q1_only = 0
    hit_but_short_q1rc_only = 0
    hit_q1rc_only = 0

    hit_but_short_q2_only = 0
    hit_q2_only = 0
    hit_but_short_q2rc_only = 0
    hit_q2rc_only = 0
    is_contam = 0
    reject_reads_list = []
    #   with pysam.FastxFile(fq_filename) as fin, gzip.open(fqout_filename + '.gz', mode='wb') as fout:
    with pysam.FastxFile(fq_filename) as fin, open(fqout_filename, mode='wt') as fout:
        for entry in fin:
            count += 1
            if decontam_tranposon == True:
                matchescontam1 = find_near_matches(contam1, entry.sequence, max_substitutions=subs,
                                                   max_deletions=dels,
                                                   max_insertions=insrt)
                matchescontam2 = find_near_matches(contam2, entry.sequence, max_substitutions=subs,
                                                   max_deletions=dels,
                                                   max_insertions=insrt)
                matchescontam1rc = find_near_matches(contam1rc, entry.sequence, max_substitutions=subs,
                                                     max_deletions=dels,
                                                     max_insertions=insrt)
                matchescontam2rc = find_near_matches(contam2rc, entry.sequence, max_substitutions=subs,
                                                     max_deletions=dels,
                                                     max_insertions=insrt)

                if matchescontam1:
                    is_contam += 1
                    reject_reads_list.append(entry.name)
                    continue

                if matchescontam2:
                    is_contam += 1
                    reject_reads_list.append(entry.name)
                    continue

                if matchescontam1rc:
                    is_contam += 1
                    reject_reads_list.append(entry.name)
                    continue

                if matchescontam2rc:
                    is_contam += 1
                    reject_reads_list.append(entry.name)
                    continue

            if fuzzy_levenshtein == False:
                matchesq1 = find_near_matches(qry1, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                              max_insertions=insrt)
                matchesq2 = find_near_matches(qry2, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                              max_insertions=insrt)
                matchesq1rc = find_near_matches(qry1rc, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                                max_insertions=insrt)
                matchesq2rc = find_near_matches(qry2rc, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                                max_insertions=insrt)
            else:
                matchesq1 = find_near_matches(qry1, entry.sequence, max_l_dist=l_dist)
                matchesq2 = find_near_matches(qry2, entry.sequence, max_l_dist=l_dist)
                matchesq1rc = find_near_matches(qry1rc, entry.sequence, max_l_dist=l_dist)
                matchesq2rc = find_near_matches(qry2rc, entry.sequence, max_l_dist=l_dist)

            # skip fastq entry if multiple matches to same motif query seq
            if (len(matchesq1) > 1):
                countqmulti += 1
                reject_reads_list.append(entry.name)
                continue
            if (len(matchesq2) > 1):
                countqmulti += 1
                reject_reads_list.append(entry.name)
                continue
            if (len(matchesq1rc) > 1):
                countqmulti += 1
                reject_reads_list.append(entry.name)
                continue
            if (len(matchesq2rc) > 1):
                countqmulti += 1
                reject_reads_list.append(entry.name)
                continue
            # skip fastq entry if multiple matches to same motif query direct and reverse complement

            if ((len(matchesq1) == 1) and (len(matchesq1rc) == 1)):
                countqqrc += 1
                reject_reads_list.append(entry.name)
                continue
            if ((len(matchesq2) == 1) and (len(matchesq2rc) == 1)):
                countqqrc += 1
                reject_reads_list.append(entry.name)
                continue
            # or matches to two incompatible motifs
            if ((len(matchesq1) == 1) and (len(matchesq2rc) == 1)):
                countqqrc += 1
                reject_reads_list.append(entry.name)
                continue
            if ((len(matchesq2) == 1) and (len(matchesq1rc) == 1)):
                countqqrc += 1
                reject_reads_list.append(entry.name)
                continue
            # process motif match pairs to extract target sequences
            if ((len(matchesq1) == 1) and (len(matchesq2) == 1)):
                countq1q2 += 1
                captured_seqstring = str(entry.sequence)[matchesq1[0].end:matchesq2[0].start]
                captured_qualstring = str(entry.quality)[matchesq1[0].end:matchesq2[0].start]
                if (matchesq2[0].start <= matchesq1[0].end):
                    wrongq2q1 += 1
                    reject_reads_list.append(entry.name)
                    continue

                if (len(captured_seqstring) >= min_length):
                    hit_q1_q2 += 1
                    fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write(captured_seqstring[0:max_length] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[0:max_length] + '\n')
                    continue
                else:
                    hit_but_short_q1_q2 += 1
                    reject_reads_list.append(entry.name)
                continue
            #            break
            if ((len(matchesq1rc) == 1) and (len(matchesq2rc) == 1)):
                countq2rcq1rc += 1
                captured_seqstring = str(entry.sequence)[matchesq2rc[0].end:matchesq1rc[0].start]
                captured_qualstring = str(entry.quality)[matchesq2rc[0].end:matchesq1rc[0].start]
                if (matchesq1rc[0].start <= matchesq2rc[0].end):
                    wrongq1rcq2rc += 1
                if (len(captured_seqstring) >= min_length):
                    hit_q2rc_q1rc += 1
                    fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write(captured_seqstring[0:max_length] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[0:max_length] + '\n')
                    continue
                else:
                    hit_but_short_q2rc_q1rc += 1
                    reject_reads_list.append(entry.name)
                continue
            # process single motif matches to extract target sequences
            if (len(matchesq1) == 1):
                countq1 += 1
                captured_seqstring = str(entry.sequence)[
                                     matchesq1[0].end:]  # nothing after colon indicates end of string
                captured_qualstring = str(entry.quality)[
                                      matchesq1[0].end:]
                if (len(captured_seqstring) >= min_length):
                    hit_q1_only += 1
                    fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write(captured_seqstring[0:max_length] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[0:max_length] + '\n')
                    continue
                else:
                    hit_but_short_q1_only += 1
                    reject_reads_list.append(entry.name)
                continue
            if (len(matchesq2rc) == 1):
                countq2rc += 1
                captured_seqstring = str(entry.sequence)[
                                     matchesq2rc[0].end:]  # nothing after colon indicates end of string
                captured_qualstring = str(entry.quality)[
                                      matchesq2rc[0].end:]
                if (len(captured_seqstring) >= min_length):
                    hit_q2rc_only += 1
                    fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write(captured_seqstring[0:max_length] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[0:max_length] + '\n')
                    continue
                else:
                    hit_but_short_q2rc_only += 1
                    reject_reads_list.append(entry.name)
                continue
            if (len(matchesq1rc) == 1):
                countq1rc += 1
                captured_seqstring = str(entry.sequence)[
                                     0:matchesq1rc[0].start]  # nothing after colon indicates end of string
                captured_qualstring = str(entry.quality)[
                                      0:matchesq1rc[0].start]
                if (len(captured_seqstring) >= min_length):
                    hit_q1rc_only += 1
                    fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write(captured_seqstring[-max_length:].translate(trans)[::-1] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[-max_length:][::-1] + '\n')
                    continue
                else:
                    hit_but_short_q1rc_only += 1
                    reject_reads_list.append(entry.name)
                continue
            if (len(matchesq2) == 1):
                countq2 += 1
                captured_seqstring = str(entry.sequence)[
                                     0:matchesq2[0].start]  # nothing after colon indicates end of string
                captured_qualstring = str(entry.quality)[
                                      0:matchesq2[0].start]
                if (len(captured_seqstring) >= min_length):
                    hit_q2_only += 1
                    fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write(captured_seqstring[-max_length:].translate(trans)[::-1] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[-max_length:][::-1] + '\n')
                    continue
                else:
                    hit_but_short_q2_only += 1
                    reject_reads_list.append(entry.name)
                continue

    # very cryptic logging needs reorganising and fixing to work with multiprocessing
    print("readcount\t", count)
    print("q1\t", countq1)
    print("q2\t", countq2)
    print("q1rc\t", countq1rc)
    print("q2rc\t", countq2rc)
    print('')
    print("q1q2\t", countq1q2)
    print("hit q1q2\t", hit_q1_q2)
    print("hit_but_short_q1_q2\t", hit_but_short_q1_q2)
    # print("q1rcq2rc\t", countq1rcq2rc)
    print('')
    print("q2rcq1rc\t", countq2rcq1rc)
    print("hit_q2rc_q1rc\t", hit_q2rc_q1rc)
    print("hit_but_short_q2rc_q1rc\t", hit_but_short_q2rc_q1rc)
    print('')
    print('')
    print("wrongq2q1\t", wrongq2q1)
    print("wrongq1rcq2rc\t", wrongq1rcq2rc)
    print("q + qrc\t", countqqrc)
    print("q|qrc multi\t", countqmulti)
    # print("contains contaminating ISS sequence tags \t", is_contam)


# def survey_fastq(resultx_reads_list, resultx_reads_dict, fqout):
def survey_fastq(resultx_reads_list, fqout):
    with pysam.FastxFile(fqout) as fh:
        for entry in fh:
            resultx_reads_list.append(entry.name)
            # resultx_reads_dict[entry.name] = len(str(entry.sequence))


flanking_fastq_result_list = []

if not nano:  # nano == False
    print(datetime.datetime.now())
    pi = multiprocessing.Pool(ncpus)
    for fq in glob.glob(fastq_dir + "*q.gz") + glob.glob(fastq_dir + "*.fastq"):
        pi.apply_async(pimms_fastq,
                       args=(fq,
                             os.path.join(out_dir,
                                          os.path.splitext(os.path.splitext(os.path.basename(fq))[0])[
                                              0] + fq_result_suffix
                                          )  # output illumina fastq output filepath
                             # (os.path.splitext(os.path.splitext(fq)[0])[0] + fq_result_suffix)))
                             )
                       )

    pi.close()
    pi.join()
    print("illumina initial PIMMS filtering completed...\n")
    print(datetime.datetime.now())

    ## match fwdrev match substrings e.g: _R1_/_R2_ --fwdrev parameter
    fqp_results_fwd = sorted(glob.glob(os.path.join(out_dir, "*" + fwdrev_wc[0] + "*" + fq_result_suffix)))
    fqp_results_rev = sorted(glob.glob(os.path.join(out_dir, "*" + fwdrev_wc[1] + "*" + fq_result_suffix)))
    print(fqp_results_fwd)
    print(fqp_results_rev)
    #

    for fwd_fqp_result, rev_fqp_result in zip(fqp_results_fwd, fqp_results_rev):
        result1_reads_list = []
        result2_reads_list = []
        print(fwd_fqp_result)
        print(rev_fqp_result)
        survey_fastq(result1_reads_list, fwd_fqp_result)
        survey_fastq(result2_reads_list, rev_fqp_result)
        a = set(result1_reads_list)
        b = set(result2_reads_list)
        c = b.difference(a)

        mrg_fqp_result_filename = re.sub(fwdrev_wc[0], '_RX_', fwd_fqp_result, count=1)  # replace  fwd substring '_R1_'
        flanking_fastq_result_list = flanking_fastq_result_list + [mrg_fqp_result_filename]
        print(mrg_fqp_result_filename)
        print(str(len(a)))
        print(str(len(b)))
        print(str(len(c)))

        with pysam.FastxFile(fwd_fqp_result) as fin, gzip.open(mrg_fqp_result_filename, mode='wt') as fout:
            for entry in fin:
                fout.write((str(entry) + '\n'))

        with pysam.FastxFile(rev_fqp_result) as fin, gzip.open(mrg_fqp_result_filename, mode='at') as fout:
            for entry in fin:
                if entry.name in c:
                    fout.write((str(entry) + '\n'))
                    # fout.write(str(entry) + '\n')

    # remove intermediate fastq files
    if parsed_args[0].rmfiles:
        deleteFileList(fqp_results_fwd)
        deleteFileList(fqp_results_rev)

    print("illumina merge of fwd/reverse data  completed...\n")


elif nano:  # nano == True
    print(datetime.datetime.now())
    pn = multiprocessing.Pool(ncpus)
    for fq in glob.glob(fastq_dir + "*q.gz") + glob.glob(fastq_dir + "*.fastq"):
        pn.apply_async(pimms_fastq,
                       args=(fq,  # input nanopore fastq filepath
                             os.path.join(out_dir,
                                          os.path.splitext(os.path.splitext(os.path.basename(fq))[0])[
                                              0] + fq_result_suffix
                                          )  # output nanopore fastq output filepath
                             )
                       #   (os.path.splitext(os.path.splitext(fq)[0])[0] + fq_result_suffix))
                       )

    pn.close()
    pn.join()
    print("nanopore PIMMS filtering completed...\n")
    print(datetime.datetime.now())

if parsed_args[0].nomap:
    print("Skipping mapping step...\n")
else:
    concat_result_fastq = concat_fastq(flanking_fastq_result_list, parsed_args[0].label[0], fq_result_suffix)
    sam_output = os.path.basename(parsed_args[0].fasta[0]) + '_' + parsed_args[0].label[0] + sam_result_suffix
    run_minimap2(concat_result_fastq, sam_output, parsed_args[0].fasta[0])

# flanking_fastq_result_list


#    result1_reads_list = []
#    result1_reads_dict = {}
# if nano == False:
# p1 = multiprocessing.Process(target=pimms_fastq, args=(fq1_filename, fqout1_filename,))
# p2 = multiprocessing.Process(target=pimms_fastq, args=(fq2_filename, fqout2_filename,))
#
# # starting process 1
# p1.start()
# # starting process 2
# p2.start()
#
# # wait until process 1 is finished
# p1.join()
# # wait until process 2 is finished
# p2.join()
#
# print(datetime.datetime.now())
#
# result1_reads_list = []
# result2_reads_list = []
# result1_reads_dict = {}
# result2_reads_dict = {}
#
# survey_fastq(result1_reads_list, result1_reads_dict, fqout1_filename)
# survey_fastq(result2_reads_list, result2_reads_dict, fqout2_filename)
#
# # ## deduplicate the resulting fq files and merge into one
# a = set(result1_reads_list)
# b = set(result2_reads_list)
# c = b.difference(a)
# u = a.union(b)
# # print(str(len(a)))
# # print(str(len(b)))
# # print(str(len(c)))
# # print(str(len(u)))
# #
#
#
# with pysam.FastxFile(fqout2_filename) as fin, open(fqoutm_filename, mode='w') as fout:
#     for entry in fin:
#         if entry.name in c:
#             fout.write(str(entry) + '\n')
#
# with pysam.FastxFile(fqout1_filename) as fin, open(fqoutm_filename, mode='a') as fout:
#     for entry in fin:
#         fout.write(str(entry) + '\n')
#
# print(str(datetime.datetime.now()) + "\tPIMMS insertion flanking matches\tfq1:" + str(len(a)) + "\tfq2:" + str(
#     len(b)) +
#       "\tcombined(non redundant):" + str(len(u)) + "\tnon overlapping:" + str(len(c)) + "\n")

# exit(0)
# fqout_stem = "PIMMS2_Test"
# fq1_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/test.IN.R1.fastq")
# fq2_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/test.IN.R2.fastq")
# fq1_filename = os.path.expanduser(
#     "~/Data/PIMMS_redo/PIMMS2_stuff/Short_read_test_data_for_2.0/PIMMS Data/PIMMS_Test_R1_001.fastp.fastq.gz")
# fq2_filename = os.path.expanduser(
#     "~/Data/PIMMS_redo/PIMMS2_stuff/Short_read_test_data_for_2.0/PIMMS Data/PIMMS_Test_R2_001.fastp.fastq.gz")
# fq1_filename = os.path.expanduser("FAL74897_pass_barcode07_d9afbb31_0.fastq.gz")
# FAL74897_pass_barcode07_d9afbb31_0.fastq.gz
# FAL74897_pass_barcode08_d9afbb31_0.fastq.gz
# FAL74897_pass_barcode08_d9afbb31_all.fastq.gz
# FAL74897_pass_barcode07_d9afbb31_all.fastq.gz
# fqout_stem = "test.IN.PIMMS_rmIS"
# fastq_dir = os.path.expanduser(
#    "~/Data/PIMMS_redo/PIMMS2_DEMO_DATA_JAN2020/Long_read_test_data_for_2.0/Native_PIMMS/UK15_Media_Input/barcode08/")
# use of os.path.join here checks and if needed adds a final path separator "/'
# fastq_dir = os.path.join(os.path.expanduser(
#    "~/Data/PIMMS_redo/PIMMS2_DEMO_DATA_JAN2020/Short_read_test_data_for_2.0/PIMMS_Data/UK15_Media_Input"), "")
# fastq_dir = "~/Data/PIMMS_redo/PIMMS2_DEMO_DATA_JAN2020/Long_read_test_data_for_2.0/Native_PIMMS/UK15_Media_Input/"
# fastq_dir = os.path.join("/Users/svzaw/Data/PIMMS_redo/PIMMS2_DEMO_DATA_JAN2020/Short_read_test_data_for_2.0/PIMMS_Data/UK15_Media_Input", '')
# fastq_dir = os.path.join(
#     "/Users/svzaw/Data/PIMMS_redo/PIMMS2_DEMO_DATA_JAN2020/Short_read_test_data_for_2.0/PIMMS_Data/UK15_Blood_Output",
#     '')


# fq1_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/Short_read_test_data_for_2.0/PIMMS Data/PIMMS_Test_R1_001.fastq.gz")
# fq2_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/Short_read_test_data_for_2.0/PIMMS Data/PIMMS_Test_R2_001.fastq.gz")
# fqout1_filename = fqout_stem + ".R1.sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".fastq"
# fqout2_filename = fqout_stem + ".R2.sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".fastq"
# fqoutm_filename = fqout_stem + ".RX.sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".fastq"
# fqout_filename = "PIMMS_Test_1_sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".R1.fastq"

# mergedregex = re.compile('(' + qry1 + ')|(' + qry2 + ')|(' + qry1rc + ')|(' + qry2rc + ')')

# print(qry1)
# print(qry1rc)
# print('')
# print(qry2)
# print(qry2rc)
# print('')
# print(str(
#     datetime.datetime.now()) + "\tFinding PIMMS insertion flanking matches...\nfq1:\t" + fq1_filename + "\nfq2:\t" + fq2_filename + "\n")
###

# count = 0
# countq1 = 0
# countq2 = 0
# countq1q2 = 0
# countq1rc = 0
# countq2rc = 0
# countq1rcq2rc = 0
# hit_but_short_q1_q2 = 0
# hit_q1_q2 = 0
# countq2rcq1rc = 0
# hit_but_short_q2rc_q1rc = 0
# hit_q2rc_q1rc = 0
# wrongq2q1 = 0
# wrongq1rcq2rc = 0
# countqqrc = 0
# countqmulti = 0
# hit_but_short_q1_only = 0
# hit_q1_only = 0
# hit_but_short_q1rc_only = 0
# hit_q1rc_only = 0
#
# hit_but_short_q2_only = 0
# hit_q2_only = 0
# hit_but_short_q2rc_only = 0
# hit_q2rc_only = 0
# ray.init()

# @ray.remote
# def rrfind_near_matches(query, entry_sequence):
#     return find_near_matches(query, entry_sequence, max_substitutions=0, max_deletions=0, max_insertions=0)
# @ray.remote
# pimms_fastq(fq1_filename, fqout1_filename)
# pimms_fastq(fq2_filename, fqout2_filename)
# print(datetime.datetime.now())
# with open('reject_reads_list.txt', 'w') as filehandle:
#     for listitem in reject_reads_list:
#         filehandle.write('%s\n' % listitem)
