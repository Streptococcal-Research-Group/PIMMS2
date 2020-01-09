from fuzzysearch import find_near_matches
import pysam
import os

trans = str.maketrans('ATGCN', 'TACGN')

qry1 = "TCAGAAAACTTTGCAACAGAACC"
qry2 = "GGTTCTGTTGCAAAGTTTAAAAA"
subs = 0
insrt = 0
dels = 0
min_length = 20
max_length = 50
max_length_index = max_length - 1
fq_filename = os.path.expanduser("~/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/test.IN.R2.fastq.gz")
fqout_filename = "test.INPY_sub" + str(subs) + "_min" + str(min_length) + "_max" + str(max_length) + ".R2.fastq"
qry1rc = qry1.translate(trans)[::-1]
qry2rc = qry2.translate(trans)[::-1]
print(qry1)
print(qry1rc)
print('')
print(qry2)
print(qry2rc)
print('')

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

with pysam.FastxFile(fq_filename) as fin, open(fqout_filename, mode='w') as fout:
    for entry in fin:
        count += 1
        matchesq1 = find_near_matches(qry1, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                      max_insertions=insrt)
        matchesq2 = find_near_matches(qry2, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                      max_insertions=insrt)
        matchesq1rc = find_near_matches(qry1rc, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                        max_insertions=insrt)
        matchesq2rc = find_near_matches(qry2rc, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                        max_insertions=insrt)
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
        # process motif matches to extract target sequences
        if ((len(matchesq1) == 1) and (len(matchesq2) == 1)):
            countq1q2 += 1
            captured_seqstring = str(entry.sequence)[matchesq1[0].end + 1: matchesq2[0].start - 1]
            captured_qualstring = str(entry.quality)[matchesq1[0].end + 1: matchesq2[0].start - 1]
            if (matchesq2[0].start <= matchesq1[0].end):
                wrongq2q1 += 1

            if (len(captured_seqstring) >= min_length):
                hit_q1_q2 += 1
                fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                fout.write(captured_seqstring[0:max_length_index] + '\n')
                fout.write('+' + '\n')
                fout.write(captured_qualstring[0:max_length_index] + '\n')
                continue
            else:
                hit_but_short_q1_q2 += 1
            continue
        #            break
        if ((len(matchesq1rc) == 1) and (len(matchesq2rc) == 1)):
            countq2rcq1rc += 1
            captured_seqstring = str(entry.sequence)[matchesq2rc[0].end + 1: matchesq1rc[0].start - 1]
            captured_qualstring = str(entry.quality)[matchesq2rc[0].end + 1: matchesq1rc[0].start - 1]
            if (matchesq1rc[0].start <= matchesq2rc[0].end):
                wrongq1rcq2rc += 1
            if (len(captured_seqstring) >= min_length):
                hit_q2rc_q1rc += 1
                fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                fout.write(captured_seqstring[0:max_length_index] + '\n')
                fout.write('+' + '\n')
                fout.write(captured_qualstring[0:max_length_index] + '\n')
                continue
            else:
                hit_but_short_q2rc_q1rc += 1
            continue
        if (len(matchesq1) == 1):
            countq1 += 1
            captured_seqstring = str(entry.sequence)[
                                 matchesq1[0].end + 1:]  # nothing after colon indicates end of string
            captured_qualstring = str(entry.quality)[
                                  matchesq1[0].end + 1:]
            if (len(captured_seqstring) >= min_length):
                hit_q1_only += 1
                fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                fout.write(captured_seqstring[0:max_length_index] + '\n')
                fout.write('+' + '\n')
                fout.write(captured_qualstring[0:max_length_index] + '\n')
                continue
            else:
                hit_but_short_q1_only += 1
            continue
        if (len(matchesq2rc) == 1):
            countq2rc += 1
            captured_seqstring = str(entry.sequence)[
                                 matchesq2rc[0].end + 1:]  # nothing after colon indicates end of string
            captured_qualstring = str(entry.quality)[
                                  matchesq2rc[0].end + 1:]
            if (len(captured_seqstring) >= min_length):
                hit_q2rc_only += 1
                fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                fout.write(captured_seqstring[0:max_length_index] + '\n')
                fout.write('+' + '\n')
                fout.write(captured_qualstring[0:max_length_index] + '\n')
                continue
            else:
                hit_but_short_q2rc_only += 1
            continue
        if (len(matchesq1rc) == 1):
            countq1rc += 1
            captured_seqstring = str(entry.sequence)[
                                 0:matchesq1rc[0].start - 1]  # nothing after colon indicates end of string
            captured_qualstring = str(entry.quality)[
                                  0:matchesq1rc[0].start - 1]
            if (len(captured_seqstring) >= min_length):
                hit_q1rc_only += 1
                fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                fout.write(captured_seqstring[-(max_length + 1):] + '\n')
                fout.write('+' + '\n')
                fout.write(captured_qualstring[-(max_length + 1):] + '\n')
                continue
            else:
                hit_but_short_q1rc_only += 1
            continue
        if (len(matchesq2) == 1):
            countq2 += 1
            captured_seqstring = str(entry.sequence)[
                                 0:matchesq2[0].start - 1]  # nothing after colon indicates end of string
            captured_qualstring = str(entry.quality)[
                                  0:matchesq2[0].start - 1]
            if (len(captured_seqstring) >= min_length):
                hit_q2_only += 1
                fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                fout.write(captured_seqstring[-(max_length + 1):] + '\n')
                fout.write('+' + '\n')
                fout.write(captured_qualstring[-(max_length + 1):] + '\n')
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
# matchesq1 = find_near_matches(qry1, trgt, max_substitutions=1, max_deletions=0, max_insertions=0)
# matchesq2 = find_near_matches(qry2, trgt, max_substitutions=1, max_deletions=0, max_insertions=0)
# matchesq1rc = find_near_matches(qry1rc, trgt, max_substitutions=1, max_deletions=0, max_insertions=0)
# matchesq2rc = find_near_matches(qry2rc, trgt, max_substitutions=1, max_deletions=0, max_insertions=0)
#
# #if (len(matches) != 1):
# #    {
# print(trgt[0:matches[0].start - 1]),
# print(trgt[matches[0].start:matches[0].end]),
# print(trgt[matches[0].end + 1:len(trgt)])
# #    }
