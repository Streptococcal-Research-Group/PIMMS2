from argparse import _SubParsersAction
import pathlib
from fuzzysearch import find_near_matches
import os
import time
import sys
import argparse
import gzip
import datetime
import multiprocessing
import glob
import pysam
from pathlib import Path


# import re
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


pimms_mls = """===========================================================================================================
Pragmatic Insertional Mutation Mapping system 2 (PIMMS2) mapping pipeline
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
python PIMMS test script.....
===========================================================================================================\n"""

# Construct the argument parser
ap = argparse.ArgumentParser(description='PIMMS2 fastq sequence processing', prog="demo_pimms2",
                             epilog="\n\n*** N.B. This is a development version ***\n \n ")
modes: _SubParsersAction = ap.add_subparsers()
# modes.required = False
findflank = modes.add_parser("find_flank", help="Mode: filter fastq files to find insertion site flanking sequence")
samcoords = modes.add_parser("sam_extract", help="Mode: extract insertion site coordinates from sam file")
otherstuff = modes.add_parser("other_stuff", help='Mode: do other good PIMMS related stuff')
# Add the arguments to the parser
findflank.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.2',
                       help="version")
findflank.add_argument("--ont", required=False, action='store_true',
                       help="nanopore reads [default: illumina paired end]")
findflank.add_argument("--single", required=False, action='store_true',
                       help="illumina single end reads (un-paired) [default: illumina paired end]")
findflank.add_argument("--rmvector", required=False, action='store_true',
                       help="attempt removal of contaminating vector")

findflank.add_argument("--lev", required=False, action='store_true',
                       help="use Levenshtein distance of 1")
findflank.add_argument("-d", "--dir", required=False, nargs=1, dest='in_dir',
                       help="directory containing input fastq files")
findflank.add_argument("-r", "--rdir", required=False, nargs=1, metavar='DIR', dest='out_dir',  # action='store',
                       default=(['pimms2_' + time.strftime("%Y%m%d_%H%M%S")]),
                       help="directory to contain fastq files results")
findflank.add_argument("--max", required=False, nargs=1, type=int,
                       help="clip results to this length [illumina:60/nano:100")
findflank.add_argument("--config", required=False, nargs=1,
                       help="read parameters from config file")
samcoords.add_argument("-s", "--sam", required=False, nargs=1,
                       help="sam file of mapped IS flanking sequences")

samcoords.add_argument("--config", required=False, nargs=1,
                       help="read parameters from config file")
otherstuff.add_argument("--config", required=False, nargs=1,
                        help="read parameters from config file")

if len(sys.argv) < 2:
    ap.print_usage()
    sys.exit(1)

# args = vars(ap.parse_args())

parsed_args = ap.parse_args()

# print((vars(parsed_args)['out_dir']))
out_dir: object = vars(parsed_args)['out_dir'][0]

if os.path.isdir(out_dir):
    print('result dir exists\n')
else:
    print('creating result dir: ' + out_dir + '\n')
    createFolder(out_dir)

# exit(0)

print(pimms_mls)

ncpus = 6
nano = False
decontam_tranposon = False
fuzzy_levenshtein = False

if decontam_tranposon == False:
    decon_tag = "nodecon"
else:
    decon_tag = "decon"

# fq_result_suffix = "_pimmsout_trim100_nodecon.fastq"

if nano:  # nano == True
    subs = 1
    l_dist = 1  # maximum Levenshtein Distance
    insrt = 1
    dels = 1
    min_length = 25
    max_length = 120
else:
    subs = 1
    l_dist = 1  # maximum Levenshtein Distance
    insrt = 0
    dels = 0
    min_length = 25
    max_length = 60

# min_length = 25
# max_length = 120

fq_result_suffix = ("_pimms2out_trim" + str(max_length) + "_" + decon_tag + ".fastq")

trans = str.maketrans('ATGCN', 'TACGN')  # complement DNA lookup

contam1 = 'GATGCTCTAGAGCATTCTCT'
contam2 = 'CATTCTCTCCATCAAGCTAT'
contam1rc = contam1.translate(trans)[::-1]  # reverse complement ([::-1] -> reverse)
contam2rc = contam2.translate(trans)[::-1]

qry1 = "TCAGAAAACTTTGCAACAGAACC"
qry2 = "GGTTCTGTTGCAAAGTTTAAAAA"
qry1rc = qry1.translate(trans)[::-1]  # reverse complement ([::-1] -> reverse)
qry2rc = qry2.translate(trans)[::-1]

max_length_index = max_length - 1
# fq1_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/test.IN.R1.fastq")
# fq2_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/test.IN.R2.fastq")
fq1_filename = os.path.expanduser(
    "~/Data/PIMMS_redo/PIMMS2_stuff/Short_read_test_data_for_2.0/PIMMS Data/PIMMS_Test_R1_001.fastp.fastq.gz")
fq2_filename = os.path.expanduser(
    "~/Data/PIMMS_redo/PIMMS2_stuff/Short_read_test_data_for_2.0/PIMMS Data/PIMMS_Test_R2_001.fastp.fastq.gz")
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
fastq_dir = os.path.join(
    "/Users/svzaw/Data/PIMMS_redo/PIMMS2_DEMO_DATA_JAN2020/Short_read_test_data_for_2.0/PIMMS_Data/UK15_Blood_Output",
    '')
# fastq_dir = os.path.join("pimms1_illpe_100k", '')

fqout_stem = "PIMMS2_Test"
# fq1_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/Short_read_test_data_for_2.0/PIMMS Data/PIMMS_Test_R1_001.fastq.gz")
# fq2_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/Short_read_test_data_for_2.0/PIMMS Data/PIMMS_Test_R2_001.fastq.gz")
fqout1_filename = fqout_stem + ".R1.sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".fastq"
fqout2_filename = fqout_stem + ".R2.sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".fastq"
fqoutm_filename = fqout_stem + ".RX.sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".fastq"
fqout_filename = "PIMMS_Test_1_sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".R1.fastq"


# mergedregex = re.compile('(' + qry1 + ')|(' + qry2 + ')|(' + qry1rc + ')|(' + qry2rc + ')')

# print(qry1)
# print(qry1rc)
# print('')
# print(qry2)
# print(qry2rc)
# print('')
# print(str(
#     datetime.datetime.now()) + "\tFinding PIMMS insertion flanking matches...\nfq1:\t" + fq1_filename + "\nfq2:\t" + fq2_filename + "\n")


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
    with pysam.FastxFile(fq_filename) as fin, open(fqout_filename, mode='w') as fout:
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

    #        print(entry.name)
    #        print(entry.sequence)
    #        print(entry.comment)
    #       print(entry.quality)
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
    print("contains contaminating ISS sequence tags \t", is_contam)


# pimms_fastq(fq1_filename, fqout1_filename)
# pimms_fastq(fq2_filename, fqout2_filename)
# print(datetime.datetime.now())
# with open('reject_reads_list.txt', 'w') as filehandle:
#     for listitem in reject_reads_list:
#         filehandle.write('%s\n' % listitem)

# def survey_fastq(resultx_reads_list, resultx_reads_dict, fqout):
def survey_fastq(resultx_reads_list, fqout):
    with pysam.FastxFile(fqout) as fh:
        for entry in fh:
            resultx_reads_list.append(entry.name)
            # resultx_reads_dict[entry.name] = len(str(entry.sequence))


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


if not nano:  # nano == False
    print(datetime.datetime.now())
    pi = multiprocessing.Pool(ncpus)
    for fq in glob.glob(fastq_dir + "*q.gz"):
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

    # fqp_results_fwd = sorted(glob.glob(fastq_dir + "*_R1_*" + fq_result_suffix))
    # fqp_results_rev = sorted(glob.glob(fastq_dir + "*_R2_*" + fq_result_suffix))
    #  fqp_results_fwd = sorted(glob.glob(out_dir + "*_R1_*" + fq_result_suffix))
    #  fqp_results_rev = sorted(glob.glob(out_dir + "*_R2_*" + fq_result_suffix))
    fqp_results_fwd = sorted(glob.glob(os.path.join(out_dir, "*_R1_*" + fq_result_suffix)))
    fqp_results_rev = sorted(glob.glob(os.path.join(out_dir, "*_R2_*" + fq_result_suffix)))
    print(fqp_results_fwd)
    print(fqp_results_rev)
    import re

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

        mrg_fqp_result_filename = re.sub('_R1_', '_RX_', fwd_fqp_result, count=1)
        print(mrg_fqp_result_filename)
        print(str(len(a)))
        print(str(len(b)))
        print(str(len(c)))

        with pysam.FastxFile(fwd_fqp_result) as fin, open(mrg_fqp_result_filename, mode='w') as fout:
            for entry in fin:
                fout.write(str(entry) + '\n')

        with pysam.FastxFile(rev_fqp_result) as fin, open(mrg_fqp_result_filename, mode='a') as fout:
            for entry in fin:
                if entry.name in c:
                    fout.write(str(entry) + '\n')


elif nano:  # nano == True
    print(datetime.datetime.now())
    pn = multiprocessing.Pool(ncpus)
    for fq in glob.glob(fastq_dir + "*q.gz"):
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

#    result1_reads_list = []
#    result1_reads_dict = {}
