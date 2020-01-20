from fuzzysearch import find_near_matches
import os
import datetime
import multiprocessing
import pysam
from pathlib import Path
# import re


# import ray

from threading import Thread

pimms_mls = """===========================================================================================================
Pragmatic Insertional Mutation Mapping system (PIMMS) mapping pipeline
===========================================================================================================
       o         o
       o         o
      //        //
     //        //
   |_||_|    |_||_|   @@@@@  @@@@@@  @@     @@  @@     @@   @@@@@@
   |@||@|    |@||@|   @@  @@   @@    @@@@ @@@@  @@@@ @@@@  @@
   |@||@|    |@||@|   @@@@@    @@    @@ @@@ @@  @@ @@@ @@   @@@
   |@||@|    |@||@|   @@       @@    @@  @  @@  @@  @  @@     @@@
   |@@@@|    |@@@@|   @@       @@    @@     @@  @@     @@       @@
   |@@@@|    |@@@@|   @@     @@@@@@  @@     @@  @@     @@  @@@@@@@
===========================================================================================================
python PIMMS test script.....
===========================================================================================================\n"""

print(pimms_mls)

trans = str.maketrans('ATGCN', 'TACGN')  # complement DNA lookup

contam1 = 'GATGCTCTAGAGCATTCTCT'
contam2 = 'CATTCTCTCCATCAAGCTAT'
contam1rc = contam1.translate(trans)[::-1]  # reverse complement ([::-1] -> reverse)
contam2rc = contam2.translate(trans)[::-1]

qry1 = "TCAGAAAACTTTGCAACAGAACC"
qry2 = "GGTTCTGTTGCAAAGTTTAAAAA"
qry1rc = qry1.translate(trans)[::-1]  # reverse complement ([::-1] -> reverse)
qry2rc = qry2.translate(trans)[::-1]
subs = 0
insrt = 0
dels = 0
min_length = 25
max_length = 50
max_length_index = max_length - 1
fq1_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/test.IN.R1.fastq")
fq2_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/test.IN.R2.fastq")
fqout_stem = "test.IN.PIMMS_rmIS"
# fqout_stem = "PIMMS_Test"
# fq1_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/Short_read_test_data_for_2.0/PIMMS Data/PIMMS_Test_R1_001.fastq.gz")
# fq2_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/Short_read_test_data_for_2.0/PIMMS Data/PIMMS_Test_R2_001.fastq.gz")
fqout1_filename = fqout_stem + ".R1.sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".fastq"
fqout2_filename = fqout_stem + ".R2.sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".fastq"
fqoutm_filename = fqout_stem + ".RX.sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".fastq"
# fqout_filename = "PIMMS_Test_1_sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".R1.fastq"


# mergedregex = re.compile('(' + qry1 + ')|(' + qry2 + ')|(' + qry1rc + ')|(' + qry2rc + ')')

print(qry1)
print(qry1rc)
print('')
print(qry2)
print(qry2rc)
print('')
print(str(
    datetime.datetime.now()) + "\tFinding PIMMS insertion flanking matches...\nfq1:\t" + fq1_filename + "\nfq2:\t" + fq2_filename + "\n")

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

    with pysam.FastxFile(fq_filename) as fin, open(fqout_filename, mode='w') as fout:
        for entry in fin:
            count += 1
            matchescontam1 = find_near_matches(contam1, entry.sequence, max_substitutions=0, max_deletions=0,
                                               max_insertions=0)
            matchescontam2 = find_near_matches(contam2, entry.sequence, max_substitutions=0, max_deletions=0,
                                               max_insertions=0)
            matchescontam1rc = find_near_matches(contam1rc, entry.sequence, max_substitutions=0, max_deletions=0,
                                                 max_insertions=0)
            matchescontam2rc = find_near_matches(contam2rc, entry.sequence, max_substitutions=0, max_deletions=0,
                                                 max_insertions=0)

            if matchescontam1:
                is_contam += 1
                continue

            if matchescontam2:
                is_contam += 1
                continue

            if matchescontam1rc:
                is_contam += 1
                continue

            if matchescontam2rc:
                is_contam += 1
                continue



            matchesq1 = find_near_matches(qry1, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                          max_insertions=insrt)
            matchesq2 = find_near_matches(qry2, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                          max_insertions=insrt)
            matchesq1rc = find_near_matches(qry1rc, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                            max_insertions=insrt)
            matchesq2rc = find_near_matches(qry2rc, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                            max_insertions=insrt)


            # matchesq1 = ray.get(rrfind_near_matches.remote(qry1, entry.sequence))
            # #print(matchesqxx)
            # matchesq2 = ray.get(rrfind_near_matches.remote(qry2, entry.sequence))
            # matchesq1rc = ray.get(rrfind_near_matches.remote(qry1rc, entry.sequence))
            # matchesq2rc = ray.get(rrfind_near_matches.remote(qry2rc, entry.sequence))
            # skip fastq entry if multiple matches to same motif query seq
            if (len(matchesq1) > 1):
                countqmulti += 1
                continue
            if (len(matchesq2) > 1):
                countqmulti += 1
                continue
            if (len(matchesq1rc) > 1):
                countqmulti += 1
                continue
            if (len(matchesq2rc) > 1):
                countqmulti += 1
                continue
            # skip fastq entry if multiple matches to same motif query direct and reverse complement

            if ((len(matchesq1) == 1) and (len(matchesq1rc) == 1)):
                countqqrc += 1
                continue
            if ((len(matchesq2) == 1) and (len(matchesq2rc) == 1)):
                countqqrc += 1
                continue
            # or matches to two incompatible motifs
            if ((len(matchesq1) == 1) and (len(matchesq2rc) == 1)):
                countqqrc += 1
                continue
            if ((len(matchesq2) == 1) and (len(matchesq1rc) == 1)):
                countqqrc += 1
                continue
            # process motif match pairs to extract target sequences
            if ((len(matchesq1) == 1) and (len(matchesq2) == 1)):
                countq1q2 += 1
                captured_seqstring = str(entry.sequence)[matchesq1[0].end:matchesq2[0].start]
                captured_qualstring = str(entry.quality)[matchesq1[0].end:matchesq2[0].start]
                if (matchesq2[0].start <= matchesq1[0].end):
                    wrongq2q1 += 1

                if (len(captured_seqstring) >= min_length):
                    hit_q1_q2 += 1
                    fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write(captured_seqstring[0:max_length] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[0:max_length] + '\n')
                    continue
                else:
                    hit_but_short_q1_q2 += 1
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


# pimms_fastq(fq1_filename, fqout1_filename)
# pimms_fastq(fq2_filename, fqout2_filename)
print(datetime.datetime.now())
p1 = multiprocessing.Process(target=pimms_fastq, args=(fq1_filename, fqout1_filename,))
p2 = multiprocessing.Process(target=pimms_fastq, args=(fq2_filename, fqout2_filename,))

# starting process 1
p1.start()
# starting process 2
p2.start()

# wait until process 1 is finished
p1.join()
# wait until process 2 is finished
p2.join()

print(datetime.datetime.now())

result1_reads_list = []
result2_reads_list = []
result1_reads_dict = {}
result2_reads_dict = {}


def survey_fastq(resultx_reads_list, resultx_reads_dict, fqout):
    with pysam.FastxFile(fqout) as fh:
        for entry in fh:
            resultx_reads_list.append(entry.name)
            resultx_reads_dict[entry.name] = len(str(entry.sequence))


survey_fastq(result1_reads_list, result1_reads_dict, fqout1_filename)
survey_fastq(result2_reads_list, result2_reads_dict, fqout2_filename)

# ## deduplicate the resulting fq files and merge into one
a = set(result1_reads_list)
b = set(result2_reads_list)
c = b.difference(a)
u = a.union(b)
# print(str(len(a)))
# print(str(len(b)))
# print(str(len(c)))
# print(str(len(u)))
#




# with pysam.FastxFile(fqout2_filename) as fin, open(fqoutm_filename, mode='w') as fout:
#     for entry in fin:
#         if entry.name in c:
#             fout.write(str(entry))


with pysam.FastxFile(fqout2_filename) as fin, open(fqoutm_filename, mode='w') as fout:
    for entry in fin:
        if entry.name in c:
            fout.write(str(entry) + '\n')

with pysam.FastxFile(fqout1_filename) as fin, open(fqoutm_filename, mode='a') as fout:
    for entry in fin:
        fout.write(str(entry) + '\n')

print(str(datetime.datetime.now()) + "\tPIMMS insertion flanking matches\tfq1:" + str(len(a)) + "\tfq2:" + str(len(b)) +
      "\tcombined(non redundant):" + str(len(u)) + "\tnon overlapping:" + str(len(c)) + "\n")
