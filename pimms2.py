import datetime
import fileinput
import glob
import gzip
import multiprocessing
import os
import random  # for log file names
import re
import shutil
import subprocess
import sys
import time
import urllib as ul  # for removing url style encoding from gff text notes
from pathlib import Path
import configargparse
import pandas as pd
import pandasql as ps

import pysam  # sequence format specific module fastq/bam/sam...
import gffpandas.gffpandas as gffpd  # annotation format specific module gff3
from fuzzysearch import find_near_matches

pimms_mssg = """
===========================================================================================================
Pragmatic Insertional Mutation Mapping system (PIMMS) mapping pipeline v2
===========================================================================================================
       o         o
       o         o
      //        //
     //        //
   |_||_|    |_||_|   @@@@@  @@@@@@  @@     @@  @@     @@   @@@@@@     @@@@@@
   |@||@|    |@||@|   @@  @@   @@    @@@@ @@@@  @@@@ @@@@  @@         @@    @@
   |@||@|    |@||@|   @@@@@    @@    @@ @@@ @@  @@ @@@ @@   @@@            @@
   |@||@|    |@||@|   @@       @@    @@  @  @@  @@  @  @@     @@@         @@
   |@@@@|    |@@@@|   @@       @@    @@     @@  @@     @@       @@      @@
   |@@@@|    |@@@@|   @@     @@@@@@  @@     @@  @@     @@  @@@@@@@    @@@@@@@@
===========================================================================================================
PIMMS2 """

pimms_mssg2 = """ mode
===========================================================================================================

"""


class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end


def create_folder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def make_results_dirs_in_sam_dir(samfile_path, run_label):
    samdirname = os.path.dirname(samfile_path)
    results_dir_db = os.path.join(samdirname, run_label + '_out_dashboard')
    results_dir_info = os.path.join(samdirname, run_label + '_out_info')
    # results_dir_kept = os.path.join(samdirname, run_label + '_out_kept')
    if not os.path.exists(results_dir_db):
        try:
            os.makedirs(results_dir_db)
            os.makedirs(results_dir_info)
        except OSError:
            print("Error while creating result dirs in {samdirname}")
    else:
        results_dir_db = os.path.join(samdirname, run_label + time.strftime("_%d%m%y_%H%M%S") + '_results_dashboard')
        results_dir_info = os.path.join(samdirname, run_label + time.strftime("_%d%m%y_%H%M%S") + '_results_info')
        try:
            os.makedirs(results_dir_db)
            os.makedirs(results_dir_info)
        except OSError:
            print("Error while creating incremented result dirs in {samdirname}")

    return samdirname, results_dir_db, results_dir_info


def delete_file_list(file_list):
    for file_path in file_list:
        try:
            os.remove(file_path)
        except OSError:
            print("Error while deleting file {filePath}")


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise configargparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def prog_in_path_check(prog_to_check):
    if shutil.which(prog_to_check):
        print(prog_to_check + ' is in the path :-)')
    else:
        sys.exit('\nERROR: ' + prog_to_check +
                 ' cannot be found in the path. \nSYS.EXIT: Please check your environment and ensure ' + prog_to_check +
                 ' is installed and available before trying again.\n\n')


def concat_fastq_raw(flanking_fastq_list, label, fq_file_suffix, concat_out_dir):
    concat_fastq_result_filename = os.path.join(concat_out_dir, label + '_RX_concat' + fq_file_suffix + '.gz')
    print(concat_fastq_result_filename)
    print(" ".join(flanking_fastq_list))

    with gzip.open(concat_fastq_result_filename, "wt", compresslevel=6) as big_file:
        with fileinput.input(files=flanking_fastq_list) as inputs:
            for line in inputs:
                big_file.write(line)

    if not parsed_args[0].keep:
        print('Removing intermediate fastq flanking reads files')
        # print(flanking_fastq_list)
        delete_file_list(flanking_fastq_list)

    return concat_fastq_result_filename


############################
# FIND_FLANK FUNCTIONS:
############################

def find_read_files_with_glob(indir, wildcards):
    for suffix_wc in wildcards:
        read_files = glob.glob(indir + suffix_wc)
        if len(read_files):
            return read_files
    sys.exit("SYS EXIT: unable to find read files, check file suffixes match permissible: " + wildcards + '\n')


def merge_logs(log_path):
    log_files = glob.glob(os.path.join(log_path, "log_*txt"))

    df_from_each_log = (pd.read_table(f) for f in log_files)
    merged_logs_df = pd.concat(df_from_each_log, ignore_index=True)
    merged_logs_df = merged_logs_df.sort_values(by=['fq_filename'])
    log_sums = merged_logs_df.sum(numeric_only=True)
    log_sums['fq_filename'] = 'COMBINED'
    merged_logs_df = merged_logs_df.append(log_sums, ignore_index=True)
    merged_logs_df.to_csv(os.path.join(log_path, "..", 'result_summary.txt'), sep='\t', index=False)
    print(merged_logs_df.to_string(index=False))
    return merged_logs_df


def run_minimap2(flanking_fastq_concat_result, sam_output_result, genome_fasta):
    stream = os.popen('minimap2 --version')
    output = stream.read()
    print('calling minimap version: ' + output)
    # process = subprocess.Popen(['minimap2', '--version'],
    print(' '.join(['minimap2', '-x', 'sr', '-a',
                    '-y',  # -y adds fastq comment to sam?
                    '-o', sam_output_result, genome_fasta, flanking_fastq_concat_result,
                    '--secondary=no', '--sam-hit-only']))
    if parsed_args[0].nano:
        mm_mode = 'map-ont'
    else:
        mm_mode = 'sr'

    process = subprocess.Popen(
        ['minimap2', '-x', mm_mode,
         '-a',
         '-y',  # -y adds fastq comment to sam
         '-o', sam_output_result,
         genome_fasta,
         flanking_fastq_concat_result,
         '--secondary=no', '--sam-hit-only'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stdout.decode('utf-8'))
    print(stderr.decode('utf-8'))


def run_bwa(flanking_fastq_concat_result, sam_output_result, genome_fasta, ncpus):
    bwa_index_dir = Path(genome_fasta).stem + '_index'
    if not os.path.exists(os.path.join(bwa_index_dir, genome_fasta + '.sa')):
        print('Creating BWA index...')
        create_folder(bwa_index_dir)
        fasta_to_index = os.path.join(bwa_index_dir, Path(genome_fasta).name)
        shutil.copyfile(genome_fasta, fasta_to_index)
        process = subprocess.Popen(
            ['bwa', 'index',
             fasta_to_index],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        print(stdout.decode('utf-8'))
        print(stderr.decode('utf-8'))
    else:
        print('Using existing BWA index...')

    print(' '.join(['bwa', 'mem', genome_fasta, flanking_fastq_concat_result, sam_output_result]))
    # with open(sam_output_result, 'w') as f:
    process = subprocess.Popen(
        ['bwa', 'mem',
         '-t', str(ncpus), "-C",
         '-o', sam_output_result,
         os.path.join(bwa_index_dir, genome_fasta),
         flanking_fastq_concat_result],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stdout.decode('utf-8'))
    print(stderr.decode('utf-8'))


def py_sam_to_bam(sam_output_result):
    bam_output_result = re.sub('.sam', '.bam', sam_output_result)
    # add?? line to remove unmapped readsbased on: | samtools view -F 4 -o onlyMapped.bam ??
    # noinspection PyUnresolvedReferences
    pysam.sort('-O' 'BAM', "-o", bam_output_result, sam_output_result)
    # noinspection PyUnresolvedReferences
    pysam.index(bam_output_result)
    print('\nMapping stats (flagstat):\n')
    # noinspection PyUnresolvedReferences
    for fsline in pysam.flagstat(bam_output_result).splitlines()[:5]:
        print(fsline)

    print('\n\n')
    # if parsed_args[0].rmfiles:
    if not parsed_args[0].keep:
        delete_file_list([sam_output_result])

    return bam_output_result


def pimms_fastq(fq_filename, fqout_filename, out_dir_logs, nano):
    trans = str.maketrans('ATGCN', 'TACGN')  # complement DNA lookup
    qry1 = parsed_args[0].motif1[0].strip("'\"")
    qry2 = parsed_args[0].motif2[0].strip("'\"")
    # print(str(qry1))
    # print(str(qry2))

    # revcomp using maketrans lookup and a string reverse
    qry1rc = qry1.translate(trans)[::-1]  # reverse complement transposon motif1 ([::-1] -> reverse)
    qry2rc = qry2.translate(trans)[::-1]  # reverse complement transposon motif2
    # print(str(qry1rc))
    # print(str(qry2rc))

    # if parsed_args[0].nano:  # nano == True
    if nano:  # nano == True
        subs = 0
        insrt = 0
        dels = 0

        fuzzy_levenshtein = True
        l_dist = parsed_args[0].lev[0]  # maximum Levenshtein Distance
        min_length = 50
        max_length = 200
        qual_char = parsed_args[0].qual_char
        print('overriding with Nanopore appropriate settings: Levenshtein distance of ' + str(
            l_dist) + ' + sequence length min = ' + str(min_length) + ', max = ' + str(max_length))
    else:
        # fuzzy_levenshtein = False
        subs = parsed_args[0].sub[0]
        l_dist = parsed_args[0].lev[0]  # maximum Levenstein Distance
        fuzzy_levenshtein = bool(l_dist)
        insrt = parsed_args[0].insert[0]
        dels = parsed_args[0].deletion[0]
        min_length = parsed_args[0].min
        max_length = parsed_args[0].max
    # print('standard settings\n')

    # print("\n" + fq_filename + "  ->>\n" + fqout_filename + "#####################\n")

    count = 0
    countq1 = 0
    countq2 = 0
    countq1q2 = 0
    countq1rc = 0
    countq2rc = 0
    # countq1rcq2rc = 0
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
    # is_contam = 0
    # reject_reads_list = []
    # reject_reads_dict = dict()

    with pysam.FastxFile(fq_filename, persist=False) as fin, open(fqout_filename, mode='wt') as fout:
        print(fq_filename, '  ==>\n\t\t\t', fqout_filename, '\n')
        for entry in fin:
            count += 1

            if not fuzzy_levenshtein:
                # print('find_near_matches \n')
                matchesq1 = find_near_matches(qry1, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                              max_insertions=insrt)
                matchesq2 = find_near_matches(qry2, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                              max_insertions=insrt)
                matchesq1rc = find_near_matches(qry1rc, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                                max_insertions=insrt)
                matchesq2rc = find_near_matches(qry2rc, entry.sequence, max_substitutions=subs, max_deletions=dels,
                                                max_insertions=insrt)
            else:
                # print('find_near_matches lev\n')
                matchesq1 = find_near_matches(qry1, entry.sequence, max_l_dist=l_dist)
                matchesq2 = find_near_matches(qry2, entry.sequence, max_l_dist=l_dist)
                matchesq1rc = find_near_matches(qry1rc, entry.sequence, max_l_dist=l_dist)
                matchesq2rc = find_near_matches(qry2rc, entry.sequence, max_l_dist=l_dist)

            if not bool(matchesq1 + matchesq2 + matchesq1rc + matchesq2rc):
                # print(matchesq1 + matchesq2 + matchesq1rc + matchesq1rc)
                # reject_reads_dict.update({entry.name: 'nomatch'})
                continue
            # skip fastq entry if multiple matches to same motif query seq
            if len(matchesq1) > 1:
                countqmulti += 1
                # reject_reads_list.append(entry.name)
                # reject_reads_dict.update({entry.name: 'multi'})
                continue
            if len(matchesq2) > 1:
                countqmulti += 1
                # reject_reads_list.append(entry.name)
                # reject_reads_dict.update({entry.name: 'multi'})
                continue
            if len(matchesq1rc) > 1:
                countqmulti += 1
                # reject_reads_list.append(entry.name)
                # reject_reads_dict.update({entry.name: 'multi'})
                continue
            if len(matchesq2rc) > 1:
                countqmulti += 1
                # reject_reads_list.append(entry.name)
                # reject_reads_dict.update({entry.name: 'multi'})
                continue
            # skip fastq entry if multiple matches to same motif query direct and reverse complement

            if (len(matchesq1) == 1) and (len(matchesq1rc) == 1):
                countqqrc += 1
                # reject_reads_list.append(entry.name)
                # reject_reads_dict.update({entry.name: 'multicomp'})
                continue
            if (len(matchesq2) == 1) and (len(matchesq2rc) == 1):
                countqqrc += 1
                # reject_reads_list.append(entry.name)
                # reject_reads_dict.update({entry.name: 'multicomp'})
                continue
            # or matches to two incompatible motifs
            if (len(matchesq1) == 1) and (len(matchesq2rc) == 1):
                countqqrc += 1
                # reject_reads_list.append(entry.name)
                # reject_reads_dict.update({entry.name: 'multicomp'})
                continue
            if (len(matchesq2) == 1) and (len(matchesq1rc) == 1):
                countqqrc += 1
                # reject_reads_list.append(entry.name)
                # reject_reads_dict.update({entry.name: 'multicomp'})
                continue
            # process motif match pairs to extract target sequences
            if (len(matchesq1) == 1) and (len(matchesq2) == 1):
                countq1q2 += 1
                captured_seqstring = str(entry.sequence)[matchesq1[0].end:matchesq2[0].start]
                captured_qualstring = str(entry.quality)[matchesq1[0].end:matchesq2[0].start]
                if len(captured_qualstring) < 5:
                    captured_qualstring = qual_char * len(captured_seqstring)

                if matchesq2[0].start <= matchesq1[0].end:
                    wrongq2q1 += 1
                    # reject_reads_list.append(entry.name)
                    # reject_reads_dict.update({entry.name: 'ooorder'})
                    continue

                if len(captured_seqstring) >= min_length:
                    hit_q1_q2 += 1
                    # print('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    # fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write('@' + str(entry.name) + ' ' + 'CO:Z:' + str(
                        entry.comment) + '\n')  # make comment bam compatible
                    fout.write(captured_seqstring[0:max_length] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[0:max_length] + '\n')
                    continue
                else:
                    hit_but_short_q1_q2 += 1
                    # reject_reads_list.append(entry.name)
                    # reject_reads_dict.update({entry.name: 'short'})
                continue
            #            break
            if (len(matchesq1rc) == 1) and (len(matchesq2rc) == 1):
                countq2rcq1rc += 1
                captured_seqstring = str(entry.sequence)[matchesq2rc[0].end:matchesq1rc[0].start]
                captured_qualstring = str(entry.quality)[matchesq2rc[0].end:matchesq1rc[0].start]
                if len(captured_qualstring) < 5:
                    captured_qualstring = qual_char * len(captured_seqstring)

                if matchesq1rc[0].start <= matchesq2rc[0].end:
                    wrongq1rcq2rc += 1
                    # reject_reads_dict.update({entry.name: 'ooorder'})
                if len(captured_seqstring) >= min_length:
                    hit_q2rc_q1rc += 1
                    # fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write('@' + str(entry.name) + ' ' + 'CO:Z:' + str(
                        entry.comment) + '\n')  # make comment bam compatible
                    fout.write(captured_seqstring[0:max_length] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[0:max_length] + '\n')
                    continue
                else:
                    hit_but_short_q2rc_q1rc += 1
                    # reject_reads_list.append(entry.name)
                    # reject_reads_dict.update({entry.name: 'short'})
                continue
            # process single motif matches to extract target sequences
            if len(matchesq1) == 1:
                countq1 += 1
                captured_seqstring = str(entry.sequence)[
                                     matchesq1[0].end:]  # nothing after colon indicates end of string
                captured_qualstring = str(entry.quality)[
                                      matchesq1[0].end:]
                if len(captured_qualstring) < 5:
                    captured_qualstring = qual_char * len(captured_seqstring)

                if len(captured_seqstring) >= min_length:
                    hit_q1_only += 1
                    # fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write('@' + str(entry.name) + ' ' + 'CO:Z:' + str(
                        entry.comment) + '\n')  # make comment bam compatible
                    fout.write(captured_seqstring[0:max_length] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[0:max_length] + '\n')
                    continue
                else:
                    hit_but_short_q1_only += 1
                    # reject_reads_list.append(entry.name)
                    # reject_reads_dict.update({entry.name: 'short'})
                continue
            if len(matchesq2rc) == 1:
                countq2rc += 1
                captured_seqstring = str(entry.sequence)[
                                     matchesq2rc[0].end:]  # nothing after colon indicates end of string
                captured_qualstring = str(entry.quality)[
                                      matchesq2rc[0].end:]
                if len(captured_qualstring) < 5:
                    captured_qualstring = qual_char * len(captured_seqstring)

                if len(captured_seqstring) >= min_length:
                    hit_q2rc_only += 1
                    # fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write('@' + str(entry.name) + ' ' + 'CO:Z:' + str(
                        entry.comment) + '\n')  # make comment bam compatible
                    fout.write(captured_seqstring[0:max_length] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[0:max_length] + '\n')
                    continue
                else:
                    hit_but_short_q2rc_only += 1
                    # reject_reads_list.append(entry.name)
                    # reject_reads_dict.update({entry.name: 'short'})
                continue
            if len(matchesq1rc) == 1:
                countq1rc += 1
                captured_seqstring = str(entry.sequence)[
                                     0:matchesq1rc[0].start]  # nothing after colon indicates end of string
                captured_qualstring = str(entry.quality)[
                                      0:matchesq1rc[0].start]
                if len(captured_qualstring) < 5:
                    captured_qualstring = qual_char * len(captured_seqstring)

                if len(captured_seqstring) >= min_length:
                    hit_q1rc_only += 1
                    # fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write('@' + str(entry.name) + ' ' + 'CO:Z:' + str(
                        entry.comment) + '\n')  # make comment bam compatible
                    fout.write(captured_seqstring[-max_length:].translate(trans)[::-1] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[-max_length:][::-1] + '\n')
                    continue
                else:
                    hit_but_short_q1rc_only += 1
                    # reject_reads_list.append(entry.name)
                    # reject_reads_dict.update({entry.name: 'short'})
                continue
            if len(matchesq2) == 1:
                countq2 += 1
                captured_seqstring = str(entry.sequence)[
                                     0:matchesq2[0].start]  # nothing after colon indicates end of string
                captured_qualstring = str(entry.quality)[
                                      0:matchesq2[0].start]
                if len(captured_qualstring) < 5:
                    captured_qualstring = qual_char * len(captured_seqstring)

                if len(captured_seqstring) >= min_length:
                    hit_q2_only += 1
                    # fout.write('@' + str(entry.name) + ' ' + str(entry.comment) + '\n')
                    fout.write('@' + str(entry.name) + ' ' + 'CO:Z:' + str(
                        entry.comment) + '\n')  # make comment bam compatible
                    fout.write(captured_seqstring[-max_length:].translate(trans)[::-1] + '\n')
                    fout.write('+' + '\n')
                    fout.write(captured_qualstring[-max_length:][::-1] + '\n')
                    continue
                else:
                    hit_but_short_q2_only += 1
                    # reject_reads_list.append(entry.name)
                    # reject_reads_dict.update({entry.name: 'short'})
                continue

    # very cryptic logging needs reorganising and fixing to work with multiprocessing
    # note added random int  to get almost unique log file names need to find fix
    log_file = open(
        f'{out_dir_logs}/log_{os.getppid()}_{multiprocessing.current_process().pid}_{random.randint(1000, 9999)}.txt', 'w'
    )

    try:
        log_file.write(
            f'fq_filename\tread count:\tmultiple copies of a motif:\tmismatched motifs:\tboth motifs (fwd|revcomp):\t'
            'both motifs (fwd|revcomp) >= {min_length}:\tsingle motif  >= {min_length}:\ttotal passed\n')
        log_file.write(f'{os.path.basename(fq_filename)}\t')
        log_file.write(f'{count}\t')
        log_file.write(f'{countqmulti}\t')
        log_file.write(f'{countqqrc}\t')
        log_file.write(
            f'{hit_q1_q2 + hit_q2rc_q1rc + hit_but_short_q1_q2 + hit_but_short_q2rc_q1rc}\t')
        log_file.write(f'{hit_q1_q2 + hit_q2rc_q1rc}\t')
        log_file.write(f'{hit_q1_only + hit_q2_only + hit_q1rc_only + hit_q2rc_only}\t')
        log_file.write(f'{hit_q1_only + hit_q2_only + hit_q1rc_only + hit_q2rc_only + hit_q1_q2 + hit_q2rc_q1rc}\n')
    except Exception as e:
        # Write problems to the error file
        log_file.write(f'ERROR: {e} problem with motif matching to {fq_filename}!\n')
    finally:
        # Close the files!
        log_file.close()


# def survey_fastq(resultx_reads_list, resultx_reads_dict, fqout):
def survey_fastq(resultx_reads_list, fqout):
    with pysam.FastxFile(fqout, persist=False) as fh:
        for entry in fh:
            resultx_reads_list.append(entry.name)
            # resultx_reads_dict[entry.name] = len(str(entry.sequence))


############################
# FIND_FLANK FUNCTIONS end
############################

############################
# SAM_COORDS FUNCTIONS:
############################


def process_gff(gff_file, gff_feat_type, gff_extra, rdir):
    annotation = gffpd.read_gff3(gff_file)
    annotation = annotation.filter_feature_of_type(gff_feat_type)

    gff_stem, gff_ext = os.path.splitext(os.path.basename(gff_file))

    if gff_feat_type[0] == "pseudogene":
        annotation.to_gff3(os.path.join(rdir, gff_stem + '_pimms_features_pseudogene.gff'))
    else:
        annotation.to_gff3(os.path.join(rdir, gff_stem + '_pimms_features.gff'))

    # break 9th gff column key=value pairs down to make additional columns
    attr_to_columns = annotation.attributes_to_columns()

    if attr_to_columns.empty:
        # return empty dataframes if no rows of required type are found print('attr_to_columns is empty!')
        return attr_to_columns, attr_to_columns

    attr_to_columns = attr_to_columns.assign(
        feat_length=(attr_to_columns.end - attr_to_columns.start + 1)).dropna(axis=1,
                                                                              how='all').drop(columns=['attributes'])

    data_top = attr_to_columns.head()
    print(data_top)

    # remove RFC 3986 % encoding from product (gff3 attribute)
    # attr_to_columns = attr_to_columns.assign(product_nopc=attr_to_columns['product'].apply(ul.parse.unquote)).drop(
    #    columns=['product']).rename(columns={'product_nopc': 'product'})

    # attr_to_columns['product'] = attr_to_columns['product'].apply(ul.parse.unquote)
    if 'product' not in attr_to_columns:
        attr_to_columns['product'] = '-'
    else:
        attr_to_columns['product'] = attr_to_columns['product'].fillna('').astype(str).apply(ul.parse.unquote)  # added fix for None datatype
    # fix to skip requested extra gff annotation field if not present in GFF
    drop_gff_extra = []
    for field in gff_extra:
        if field not in attr_to_columns:
            print("Warning: Unable to find '" + field + "' in " + str(gff_file) + ' file, continuing...')
            drop_gff_extra.append(field)

    gff_extra = [item for item in gff_extra if item not in drop_gff_extra]

    # Remove URL character encoding from columns  (skipping translation if present as this breaks the decoding
    for field in gff_extra:
        if field == 'translation':
            continue
        else:
            attr_to_columns[field] = attr_to_columns[field].fillna('').astype(str).apply(ul.parse.unquote)  # added fix for None datatype

    gff_columns_addback = attr_to_columns[['seq_id',
                                           'ID',  # additional hopefully unique feature ID
                                           'locus_tag',
                                           'type',
                                           'gene',
                                           'start',
                                           'end',
                                           'feat_length',
                                           'product'] + gff_extra].copy()  # add extra fields from gff
    # fix to remove na values and allow joining with processed data also processed with fillna to allow group_by usage
    # note .copy on previous line
    gff_columns_addback.fillna('', inplace=True)

    data_top = gff_columns_addback.head()
    print(data_top)
    data_top2 = attr_to_columns.head()
    print(data_top2)

    return gff_columns_addback, attr_to_columns
    # end process_gff()


def modify_sam_stem(sam_file, min_depth_cutoff, fraction_mismatch):
    sam_stem, sam_ext = os.path.splitext(os.path.basename(sam_file))
    sam_stem: str = sam_stem + '_md' + str(min_depth_cutoff) + '_mm' + str(fraction_mismatch or '0')
    return sam_stem


def process_sam(sam_file, min_depth_cutoff, fraction_mismatch):
    sam_stem = modify_sam_stem(sam_file, min_depth_cutoff, fraction_mismatch)
    samfile = pysam.AlignmentFile(sam_file)  # without , "rb" should auto detect sam or bams

    open(sam_stem + ".bed", 'w').close()
    f = open(sam_stem + ".bed", "a")

    strand = ["+", "-"]
    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        read_str = strand[int(read.is_reverse)]
        read_bed = [read.reference_name, read.pos, read.reference_end, ".", read.mapping_quality, read_str, '# ' + read.query_name]  # read group added
        f.write('\t'.join([str(i) for i in read_bed]))
        f.write('\n')
    f.close()
    print('#BED')
    samfile.close()

    samfile = pysam.AlignmentFile(sam_file)
    open(sam_stem + "_insert_coords.txt", 'w').close()
    f2 = open(sam_stem + "_insert_coords.txt", "a")
    f2.write('\t'.join([str(i) for i in ['ref_name', 'coord', 'strand', 'read_name', 'read_grp', 'read_comment']]))
    # f2.write('\t'.join([str(i) for i in ['ref_name', 'coord', 'strand', 'read_name', 'read_grp']]))
    f2.write('\n')
    strand = ["+", "-"]
    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        # print(read.query_name + '.')
        # continue
        nm_value = read.get_tag('NM')
        if fraction_mismatch:  # and NM_value > 0:
            if (read.query_alignment_length * fraction_mismatch[0]) > nm_value:
                continue
        read_str = strand[int(read.is_reverse)]  # coverts is_reverse boolean into + or - strings
        # print(STR)
        # continue
        if read_str == '+':
            read_coords = [read.reference_name, (read.reference_start + 4), read_str, '# ' + read.query_name, ':'.join(read.query_name.split(':', 4)[:4]),
                           read.get_tag("CO").split(":")[-1]]  # add fq comment sample id number
            #          ':'.join(read.query_name.split(':', 4)[:4])]
            f2.write('\t'.join([str(i) for i in read_coords]))
            f2.write('\n')
        if read_str == '-':
            read_coords = [read.reference_name, (read.reference_end - 4), read_str, '# ' + read.query_name, ':'.join(read.query_name.split(':', 4)[:4]),
                           read.get_tag("CO").split(":")[-1]]  # add fq comment sample id number
            #          ':'.join(read.query_name.split(':', 4)[:4])]
            f2.write('\t'.join([str(i) for i in read_coords]))
            f2.write('\n')
    f2.close()
    samfile.close()
    print('#COORDS')
    # end process_sam()


def seqid_consistancy_check(mygffcolumns, my_sam):
    af = pysam.AlignmentFile(my_sam)
    sam_seq_id_list = [name['SN'] for name in af.header['SQ']]
    gff_seq_id_list = mygffcolumns.seq_id.unique().tolist()
    sam_seq_id_list.sort()
    gff_seq_id_list.sort()
    if sam_seq_id_list == gff_seq_id_list:
        print('GFF & mapping reference sequence IDs match')
    elif parsed_args[0].gff_force:
        print('\nWARNING: GFF & mapping reference sequence IDs are inconsistent. \n' +
              'sequence ID mismatch overridden by --gff_force\ngff:\n' +
              str(gff_seq_id_list) + '\nsam/bam:\n' + str(sam_seq_id_list) + '\n')
    else:
        sys.exit(
            '\nERROR: GFF & mapping reference sequence IDs are inconsistent. \n' +
            'SYS.EXIT: Please check and update the sequence IDs in your sequence and gff files so they match up before running again.\ngff:\n' +
            str(gff_seq_id_list) + '\nsam/bam:\n' + str(sam_seq_id_list) + '\n' +
            'NOTE: If the sequence ID mismatch is benign e.g. an extra plasmid/contig, override by using --gff_force with bam_extract/full_process\n')

    print(type(sam_seq_id_list))
    print(type(gff_seq_id_list))
    print(sam_seq_id_list)
    print(gff_seq_id_list)


def coordinates_to_features_reps(sam_stem, attr_to_columns, condition_label):
    coord_reps_df = pd.read_csv(sam_stem + "_insert_coords.txt", sep='\t', dtype={'ref_name': "str",
                                                                                  'coord': "int64",
                                                                                  'read_comment': "str",
                                                                                  # adding miseq support
                                                                                  'read_grp': "str"})

    read_grps = sorted(coord_reps_df.read_grp.unique())
    read_comments = sorted(coord_reps_df.read_comment.unique())
    print(str(len(read_grps)) + ' readgroups found:')
    print('\n'.join(read_grps))
    print(str(len(read_comments)) + ' sample comments found:')
    print(', '.join(read_comments))

    if (len(read_grps) < len(read_comments)) & (len(read_comments) >= 3):
        coord_counts_reps_df = coord_reps_df.drop('read_grp', 1).rename(columns={"read_comment": 'sample_info'},
                                                                        inplace=False).groupby(["ref_name",
                                                                                                "coord",
                                                                                                'sample_info']).size().reset_index(
            name=condition_label + '_')  # adding miseq support
        print(str(len(read_comments)) + " sample replicates/mutant pools established")
        # print(coord_counts_reps_df.head())

    elif (len(read_grps) >= 3) & (len(read_comments) < len(read_grps)):
        coord_counts_reps_df = coord_reps_df.drop('read_comment', 1).rename(columns={"read_grp": 'sample_info'},
                                                                            inplace=False).groupby(["ref_name",
                                                                                                    "coord",
                                                                                                    "sample_info"]).size().reset_index(
            name=condition_label + '_')  # adding miseq support
        print(str(len(read_grps)) + " sample replicates/mutant pools established")
        # print(coord_counts_reps_df.head())

    # if max(len(read_grps), len(read_comments)) < 3:
    else:
        print(
            "Warning: Unable to resolve >= 3 samples in fastq/sam/bam data, continuing without replicate insertion counts" +
            "\nN.B: If this is an error the software may need updating to recognise novel fastq naming conventions")
        # returning an empty dataframe
        return pd.DataFrame()

    coord_df_pivot = coord_counts_reps_df.copy(deep=False).pivot_table(index=["ref_name", "coord"],
                                                                       columns=['sample_info'],
                                                                       values=[condition_label + '_'],
                                                                       fill_value=0).reset_index()

    coord_df_pivot.columns = [''.join(col).strip() for col in coord_df_pivot.columns.values]

    sample_grps = sorted(coord_counts_reps_df.sample_info.unique())

    old_rep_names = [condition_label + '_' + str(x) for x in sample_grps]
    new_rep_names = [condition_label + '_' + "MP" + str(x) for x in range(1, len(sample_grps) + 1)]

    coord_df_pivot.rename(columns=dict(zip(old_rep_names, new_rep_names)), inplace=True)

    attr_to_columns_short = attr_to_columns[["seq_id", "start", "end"]]

    sqlcode = '''
        select coord_df_pivot.*
        ,attr_to_columns_short.*
        from attr_to_columns_short
        left join coord_df_pivot
        on coord_df_pivot.coord between attr_to_columns_short.start and attr_to_columns_short.end
        where coord_df_pivot.ref_name like '%' || attr_to_columns_short.seq_id || '%'
        '''
    # wierd sqlite concatenation + >> || ,  '%' == wildcard double check effect of this
    # this line should allow multi contig files

    mp_coords_join_gff = ps.sqldf(sqlcode, locals())

    # remove first 2 columns ref_name, coord sums the rest according to the feature coordinate groups
    mp_reps_feature_counts = mp_coords_join_gff.drop(mp_coords_join_gff.columns[[0, 1]], axis=1).groupby(
        ['seq_id', 'start', 'end']).agg(["sum"]).reset_index()

    mp_reps_feature_counts.columns = mp_reps_feature_counts.columns.get_level_values(0)
    return mp_reps_feature_counts


def coordinates_to_features(sam_stem, attr_to_columns, gff_columns_addback, condition_label, min_depth_cutoff,
                            gff_extra, db_rdir):
    coord_df = pd.read_csv(sam_stem + "_insert_coords.txt", sep='\t', dtype={'ref_name': "str", 'coord': "int64"})
    coord_counts_df = coord_df.groupby(['ref_name', 'coord']).size().reset_index(name='counts')
    # print(coord_counts_df.head())
    number_of_insertion_sites = len(coord_counts_df)
    number_of_reads_mapped = coord_counts_df['counts'].sum()
    min_reads_at_site = coord_counts_df['counts'].min()
    max_reads_at_site = coord_counts_df['counts'].max()
    median_reads_at_site = round(coord_counts_df['counts'].median(), 2)

    mean_insertion_site_depth = round(number_of_reads_mapped / number_of_insertion_sites, 2)

    coord_counts_df = coord_counts_df[coord_counts_df['counts'] >= min_depth_cutoff]

    # format insertion site info as a GFF
    coord_counts_df_pimms2_gff = coord_counts_df.reset_index()

    coord_counts_df_pimms2_gff['source'] = 'pimms2'
    coord_counts_df_pimms2_gff['feature_type'] = 'misc_feature'
    coord_counts_df_pimms2_gff['strand'] = '.'
    coord_counts_df_pimms2_gff['phase'] = '.'
    coord_counts_df_pimms2_gff['stop'] = coord_counts_df_pimms2_gff['coord']
    coord_counts_df_pimms2_gff = coord_counts_df_pimms2_gff.rename(columns={'counts': 'score', 'coord': 'start'})
    coord_counts_df_pimms2_gff['info'] = 'note=insertion;'
    coord_counts_df_pimms2_gff = coord_counts_df_pimms2_gff[
        ['ref_name', 'source', 'feature_type', 'start', 'stop', 'score', 'strand', 'phase', 'info']]
    print(coord_counts_df_pimms2_gff.head())
    coord_counts_df_pimms2_gff.to_csv(os.path.join(db_rdir, condition_label + "_pimms_insert_coordinates" + ".gff"), index=False, sep='\t',
                                      header=False)

    # added .loc to fix warning
    # SettingWithCopyWarning:
    # A value is trying to be set on a copy of a slice from a DataFrame.
    # Try using .loc[row_indexer,col_indexer] = value instead
    # coord_counts_df = \
    coord_counts_df.loc[:, 'between_insertion_gap'] = coord_counts_df.groupby(['ref_name'])['coord'].diff()
    # coord_counts_df = coord_counts_df.loc[:, 'between_insertion_gap'] = coord_counts_df['coord'].diff()
    # coord_counts_gt1_df = coord_counts_gt1_df({'between_insertion_gap': 0})

    min_between_insertion_gap = coord_counts_df['between_insertion_gap'].min()
    max_between_insertion_gap = coord_counts_df['between_insertion_gap'].max()
    median_between_insertion_gap = coord_counts_df['between_insertion_gap'].median()
    mean_between_insertion_gap = round(coord_counts_df['between_insertion_gap'].mean(), 2)

    sqlcode = '''
        select coord_counts_df.*
        ,attr_to_columns.*
        from attr_to_columns
        left join coord_counts_df
        on coord_counts_df.coord between attr_to_columns.start and attr_to_columns.end
        where coord_counts_df.ref_name like '%' || attr_to_columns.seq_id || '%'
        '''
    # wierd sqlite concatenation + >> || ,  '%' == wildcard double check effect of this
    # this line should allow multi contig files

    coords_join_gff = ps.sqldf(sqlcode, locals())

    # debugging save of intermediate data
    # coords_join_gff.to_csv("pimms_coords_join_gffstuff" + condition_label + ".txt", index=False, sep='\t', header=False)

    # quick dataframe summary
    # coords_join_gff.count
    # add position as percentile (needs manual confirmation)
    coords_join_gff = coords_join_gff.assign(
        # python/pandas implementation of PIMMS.pl code to derive insert position as percentile of gene length
        # sprintf("%.1f", ((($in-$in_start)+1) / ($in_length/100)));
        # sprintf("%.1f", ((($in_stop-$ in) + 1) / ($in_length / 100)));
        posn_as_percentile=(((coords_join_gff.coord - coords_join_gff.start) + 1) / (
                coords_join_gff.feat_length / 100)).where(
            coords_join_gff.strand == '+', ((coords_join_gff.end - coords_join_gff.coord) + 1) / (
                    coords_join_gff.feat_length / 100))).round({"posn_as_percentile": 1})
    print(list(attr_to_columns.columns.values))

    # Important fix group by doesn't work -- any rows with nan values get dropped *yikes very bad!!!!*
    coords_join_gff.fillna('', inplace=True)

    pimms_result_table = coords_join_gff.groupby(
        ['seq_id', 'ID',  # added ID field as unique identifier
         'locus_tag', 'type', 'gene', 'start', 'end', 'feat_length', 'product'] + gff_extra).agg(
        num_insertions_mapped_per_feat=('counts', 'sum'),
        num_insert_sites_per_feat=('counts', 'count'),
        first_insert_posn_as_percentile=('posn_as_percentile', 'min'),
        last_insert_posn_as_percentile=('posn_as_percentile', 'max')
    ).reset_index()

    # test diagnostic files
    # pimms_result_table.to_csv("pimms_coords_join_prt1_" + condition_label + ".txt", index=False, sep='\t', header=True)

    pimms_result_table = pimms_result_table.assign(num_insert_sites_per_feat_per_kb=(
            (pimms_result_table.num_insert_sites_per_feat / pimms_result_table.feat_length) * 1000),
        NRM_score=((pimms_result_table.num_insertions_mapped_per_feat / (
                pimms_result_table.feat_length / 1000)) / (
                           number_of_reads_mapped / 1e6)),
        NIM_score=((pimms_result_table.num_insert_sites_per_feat / (
                pimms_result_table.feat_length / 1000)) / (
                           number_of_reads_mapped / 1e6))
    ).round({'num_insert_sites_per_feat_per_kb': 2, 'NRM_score': 2, 'NIM_score': 2})
    print(list(pimms_result_table.columns.values))

    # test diagnostic files
    # pimms_result_table.to_csv("pimms_coords_join_prt2_" + condition_label + ".txt", index=False, sep='\t', header=False)

    pimms_result_table = pimms_result_table[['seq_id',
                                             'ID',
                                             'locus_tag',
                                             'type',
                                             'gene',
                                             'start',
                                             'end',
                                             'feat_length',
                                             'product'] + gff_extra +
                                            ['num_insertions_mapped_per_feat',
                                             'num_insert_sites_per_feat',
                                             'num_insert_sites_per_feat_per_kb',
                                             'first_insert_posn_as_percentile',
                                             'last_insert_posn_as_percentile',
                                             'NRM_score',  # Normalised Reads Mapped
                                             'NIM_score']]  # Normalised Insertions Mapped

    print(list(pimms_result_table.columns.values))

    # pimms_result_table_full gff_columns_addback

    navalues = {'num_insertions_mapped_per_feat': int(0),
                'num_insert_sites_per_feat': int(0),
                'num_insert_sites_per_feat_per_kb': int(0),
                'first_insert_posn_as_percentile': int(0),
                'last_insert_posn_as_percentile': int(0),
                'NRM_score': int(0),
                'NIM_score': int(0)}
    pimms_result_table_full = pd.merge(gff_columns_addback, pimms_result_table, how='left').fillna(value=navalues)

    # test diagnostic files
    # gff_columns_addback.to_csv("pimms_coords_join_gff_columns_addback_" + condition_label + ".txt", index=False, sep='\t', header=False)
    # pimms_result_table_full.to_csv("pimms_coords_join_prtf1_" + condition_label + ".txt", index=False, sep='\t', header=False)

    # if set add prefix to columns
    if condition_label:
        label_cols = pimms_result_table_full.columns[
            pimms_result_table_full.columns.isin(['num_insertions_mapped_per_feat',
                                                  'num_insert_sites_per_feat',
                                                  'num_insert_sites_per_feat_per_kb',
                                                  'first_insert_posn_as_percentile',
                                                  'last_insert_posn_as_percentile',
                                                  'NRM_score',
                                                  'NIM_score'])]
        pimms_result_table_full.rename(columns=dict(zip(label_cols, condition_label + '_' + label_cols)),
                                       inplace=True)

    return pimms_result_table_full
    # end of coordinates_to_features()


############################
# SAM_COORDS FUNCTIONS end
############################

def parse_arguments():
    ap = configargparse.ArgumentParser(  # description='PIMMS2 sam/bam processing',
        prog="pimms2",
        add_config_file_help=False,
        config_file_parser_class=configargparse.DefaultConfigFileParser,
        epilog="\n\n*** N.B. This is a development version ***\n \n ",
        description='''description here'''
    )
    ap.add_argument('-v', '--version', action='version', version='%(prog)s 2.0.7 demo')

    modes = ap.add_subparsers(parser_class=configargparse.ArgParser, dest='command')

    findflank = modes.add_parser("find_flank", add_config_file_help=False,
                                 help="Mode: find read regions flanking the IS sequence by mapping them to the target genome",
                                 description="Args that start with '--' (eg. --fasta) can also be set in a config file (specified via -c)")

    samcoords = modes.add_parser("bam_extract", add_config_file_help=False,
                                 help="Mode: extract insertion site coordinates from sam file",
                                 description="Args that start with '--' (eg. --fasta) can also be set in a config file (specified via -c)")

    tablemerge = modes.add_parser("table_merge", add_config_file_help=False,
                                  help='Mode: merge two compatible PIMMS results tables '
                                       '(N.B: this step does a simple table join and does not check the data)').add_mutually_exclusive_group()

    fullprocess = modes.add_parser("full_process", add_config_file_help=False,
                                   help="Mode: find_flank + bam_extract",
                                   description="Args that start with '--' (eg. --fasta) can also be set in a config file (specified via -c)")

    # FIND_FLANK args
    # to fix: nargs='?' deal with mistaken use of nargs=1 which give a single element list
    findflank.add_argument("-c", "--config", required=False, is_config_file=True,  # dest='config_file',
                           metavar='pimms2.config',
                           help="use parameters from config file")
    findflank.add_argument("--nano", required=False, action='store_true', default=False,
                           help="override with settings more suitable for nanopore")
    findflank.add_argument("--fasta", required=False, nargs=1, metavar='ref_genome.fasta', type=extant_file,
                           help="fasta file for reference genome ")
    findflank.add_argument("--qual_char", required=False, nargs='?', type=str, default='0', choices=[chr(x + 33) for x in list(range(12, 31))],
                           help="substitute a quality score ascii character when fasta read files used (nanopore only) (phred +33: ascii +:?) ['0']")
    findflank.add_argument("--nomap", required=False, action='store_true', default=False,
                           help="do not run mapping step")
    findflank.add_argument("--mapper", required=False, nargs='?', type=str, default='bwa', choices=['minimap2', 'bwa'],
                           help="select mapping software from available options")
    # findflank.add_argument("--rmfiles", required=False, action='store_true', default=False,
    #                       help="remove intermediate files")
    findflank.add_argument("--keep", required=False, action='store_true', default=False,
                           help="keep intermediate fastq files etc for diagnostic purposes")
    findflank.add_argument("--lev", required=False, nargs=1, type=int, default=0,
                           help="use Levenshtein distance (combined insert|del|sub score)")
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
    findflank.add_argument("--out_dir", required=False, nargs=1, metavar='out_dir', default=[''],
                           action='store',
                           help="directory to contain result files ['pimms2_`label`_`dmy`_`HMS`']")
    findflank.add_argument("--cpus", required=False, nargs=1, type=int,  # default=[4],
                           default=[int(os.cpu_count() / 2)],
                           help="number of processors to use [(os.cpu_count() / 2)] ")
    findflank.add_argument("--max", required=False, nargs=1, type=int, default=60,
                           help="clip results to this length [illumina:60/nano:100]")
    findflank.add_argument("--min", required=False, nargs=1, type=int, default=25,
                           help="minimum read length [illumina:60/nano:100]")
    findflank.add_argument("--motif1", required=False, nargs=1, type=str, default=['TCAGAAAACTTTGCAACAGAACC'],
                           # revcomp: GGTTCTGTTGCAAAGTTTTCTGA
                           help="IS end reference motif1 [TCAGAAAACTTTGCAACAGAACC](pGh9)")
    findflank.add_argument("--motif2", required=False, nargs=1, type=str, default=['GGTTCTGTTGCAAAGTTTAAAAA'],
                           # revcomp: TTTTTAAACTTTGCAACAGAACC
                           help="IS end reference motif2 [GGTTCTGTTGCAAAGTTTAAAAA](pGh9)")
    findflank.add_argument("--label", required=True, nargs=1, metavar='condition_name', default=[''],
                           help="identifying text tag to add to results file")

    # SAM EXTRACT args
    samcoords.add_argument("-c", "--config", required=False, is_config_file=True,
                           metavar='pimms2.config',
                           help="use parameters from config file")
    samcoords.add_argument("--bam", required=True, nargs=1, metavar='pimms.bam/sam', type=extant_file,
                           help="bam/sam file of mapped IS flanking sequences ")
    samcoords.add_argument("--nano", required=False, action='store_true', default=False,
                           help="override with settings more suitable for nanopore")
    samcoords.add_argument("--label", required=True, nargs=1, metavar='condition_name', default=[''],
                           help="text tag to add to results file")
    samcoords.add_argument("--mismatch", required=False, nargs=1, type=float, metavar='float', default=[None],
                           choices=[round(x * 0.01, 2) for x in range(0, 21)],
                           help="fraction of permitted mismatches in mapped read ( 0 <= mismatch < 0.2) [no filter]")
    samcoords.add_argument("--min_depth", required=False, nargs=1, type=int, default=2, metavar='int',
                           help="minimum read depth at insertion site >= int [2]")
    samcoords.add_argument("--noreps", required=False, action='store_true', default=False,
                           help="do not separate illumina read groups as replicate insertion count columns")
    samcoords.add_argument("--gff", required=True, nargs=1, type=extant_file, default='', metavar='genome.gff',
                           help="GFF3 formatted file to use\n(note fasta sequence present in the file must be deleted before use)")
    samcoords.add_argument("--gff_extra", required=False, nargs=1, type=str, default='', metavar="'x,y,z'",
                           help="comma separated list of extra fields to include from the GFF3 annotation\ne.g. 'ID,translation,note' ")
    samcoords.add_argument("--gff_force", required=False, action='store_true', default=False,
                           help="override GFF/BAM seq id discrepancies e.g. use when the gff has a plasmid not present in the reference sequence or vice-versa")
    samcoords.add_argument("--out_fmt", required=False, nargs=1, type=str, default=['xlsx'],
                           choices=['xlsx', 'tsv', 'csv'],
                           help="set results table file format tab/comma separated or Excel (tsv|csv|xlsx) [xlsx]")

    # TABLE_MERGE args
    tablemerge.add_argument("--xlsx", required=False, nargs=2, type=extant_file,
                            help="2x .xlsx Excel files")
    tablemerge.add_argument("--csv", required=False, nargs=2, type=extant_file,
                            help="2x .csv comma separated text/table files")
    tablemerge.add_argument("--tsv", required=False, nargs=2, type=extant_file,
                            help='2x .tsv tab (\\t) separated text/table files')

    # FULL_PROCESS ##########################
    fullprocess.add_argument("-c", "--config", required=False, is_config_file=True,
                             metavar='pimms2_run.config',
                             help="use parameters from config file")
    fullprocess.add_argument("--nano", required=False, action='store_true', default=False,
                             help="override with settings more suitable for nanopore")
    fullprocess.add_argument("--qual_char", required=False, nargs='?', type=str, default='0', choices=[chr(x + 33) for x in list(range(12, 31))],
                             help="substitute a quality score ascii character when fasta read files used (nanopore only) (phred +33: ascii +:?) ['0']")
    fullprocess.add_argument("--fasta", required=False, nargs=1, metavar='ref_genome.fasta', type=extant_file,
                             help="fasta file for reference genome ")
    fullprocess.add_argument("--nomap", required=False, action='store_true', default=False,
                             help="do not run mapping step")
    fullprocess.add_argument("--mapper", required=False, nargs='?', type=str, default='bwa', choices=['minimap2', 'bwa'],
                             help="select mapping software from available options")
    fullprocess.add_argument("--keep", required=False, action='store_true', default=False,
                             help="keep intermediate files for diagnostic purposes")
    fullprocess.add_argument("--lev", required=False, nargs=1, type=int, default=0,
                             help="use Levenshtein distance (combined insert|del|sub score)")
    fullprocess.add_argument("--sub", required=False, nargs=1, type=int, default=1,
                             help="number of permitted base substitutions in motif match [1]")
    fullprocess.add_argument("--insert", required=False, nargs=1, type=int, default=0,
                             help="number of permitted base insertions in motif match [0]")
    fullprocess.add_argument("--del", required=False, nargs=1, type=int, default=0, dest='deletion',
                             help="number of permitted base insertions in motif match [0]")
    fullprocess.add_argument("--in_dir", required=True, nargs=1, dest='in_dir', type=extant_file,
                             help="directory containing input fastq files (assumed to match '*q.gz' or '*.fastq')")
    fullprocess.add_argument("--fwdrev", required=False, nargs=1, type=str, default=['_R1_,_R2_'],
                             help="text substring to uniquely identify illumina fwd/rev paired fastq files ['_R1_,_R2_']")
    fullprocess.add_argument("--out_dir", required=False, nargs=1, metavar='out_dir', default=[''],
                             action='store',
                             help="directory to contain result files ['pimms2_`label`_`dmy`_`HMS`']")
    fullprocess.add_argument("--cpus", required=False, nargs=1, type=int,  # default=int(4),
                             default=[int(os.cpu_count() / 2)],
                             help="number of processors to use [(os.cpu_count() / 2)] ")
    fullprocess.add_argument("--max", required=True, nargs=1, type=int, default=60,
                             help="clip results to this length [illumina:60/nano:100]")
    fullprocess.add_argument("--min", required=True, nargs=1, type=int, default=25,
                             help="minimum read length [illumina:60/nano:100]")
    fullprocess.add_argument("--motif1", required=False, nargs=1, type=str, default=['TCAGAAAACTTTGCAACAGAACC'],
                             # revcomp: GGTTCTGTTGCAAAGTTTTCTGA
                             help="IS end reference motif1 [TCAGAAAACTTTGCAACAGAACC](pGh9)")
    fullprocess.add_argument("--motif2", required=False, nargs=1, type=str, default=['GGTTCTGTTGCAAAGTTTAAAAA'],
                             # revcomp: TTTTTAAACTTTGCAACAGAACC
                             help="IS end reference motif2 [GGTTCTGTTGCAAAGTTTAAAAA](pGh9)")
    fullprocess.add_argument("--label", required=True, nargs=1, metavar='condition_name', default=[''],
                             help="identifying text tag to add to results file")
    fullprocess.add_argument("--bam", required=False, nargs=1, metavar='pimms.bam/sam',  # type=extant_file,
                             type=str, default=['bam?'],
                             help=configargparse.SUPPRESS)
    # samcoords.add_argument("--nano", required=False, action='store_true', default=False,
    #                       help="override with settings more suitable for nanopore")
    # samcoords.add_argument("--label", required=False, nargs=1, metavar='condition_name', default=[''],
    #                       help="text tag to add to results file")
    fullprocess.add_argument("--mismatch", required=False, nargs=1, type=float, metavar='float', default=[None],
                             choices=[round(x * 0.01, 2) for x in range(0, 21)],
                             help="fraction of permitted mismatches in mapped read ( 0 <= mismatch < 0.2) [no filter]")
    fullprocess.add_argument("--min_depth", required=False, nargs=1, type=int, default=2, metavar='int',
                             help="minimum read depth at insertion site >= int [2]")
    fullprocess.add_argument("--noreps", required=False, action='store_true', default=False,
                             help="do not separate illumina read groups as replicate insertion count columns")
    fullprocess.add_argument("--gff", required=True, nargs=1, type=extant_file, default='', metavar='genome.gff',
                             help="GFF3 formatted file to use\n(note fasta sequence present in the file must be deleted before use)")
    fullprocess.add_argument("--gff_extra", required=False, nargs=1, type=str, default='', metavar="'x,y,z'",
                             help="comma separated list of extra fields to include from the GFF3 annotation\ne.g. 'ID,translation,note' ")
    fullprocess.add_argument("--gff_force", required=False, action='store_true', default=False,
                             help="override GFF/BAM seq id discrepancies "
                                  "e.g. use when the gff has a plasmid not present in the reference sequence or vice-versa")
    fullprocess.add_argument("--out_fmt", required=False, nargs=1, type=str, default=['xlsx'],
                             choices=['xlsx', 'tsv', 'csv'],
                             help="set results table file format tab/comma separated or Excel (tsv|csv|xlsx) [xlsx]")

    local_parsed_args = ap.parse_known_args()
    print("-----------------")
    print(local_parsed_args)
    print("-----------------")
    print(ap.format_values())
    print("-----------------")

    # exit and print short help message if no mode/arguments supplied
    if len(sys.argv) <= 2:
        ap.print_usage()
        sys.exit(1)

    if local_parsed_args[0].command == 'find_flank':
        if not local_parsed_args[0].nomap:
            prog_in_path_check(local_parsed_args[0].mapper)
            # prog_in_path_check('bwa')
            if local_parsed_args[0].fasta is None:
                ap.error("unless the --nomap flag is used please supply a sequence file e.g:  --fasta contigs.fasta")
            elif not local_parsed_args[0].label:
                ap.error("unless the --nomap flag is used please supply a text label string  --label cond_01")
            else:
                print("refseq provided: " + local_parsed_args[0].fasta[0])

    # print("##########")
    # print(ap.format_values())  # useful for logging where different settings came from
    # sys.exit(1)
    # print("\n\n\n")
    # print(parsed_args[0].command)
    # print("----------======")
    # print(ap.)
    # sys.exit(1)

    return local_parsed_args


# elif parsed_args[0].command == 'bam_extract':
def bam_extract_func(parsed_args_be):
    print(pimms_mssg + parsed_args_be[0].command + pimms_mssg2)

    if parsed_args_be[0].nano:
        parsed_args_be[0].noreps = True
    # sort out extra requested gff annotation fields
    if parsed_args_be[0].gff_extra:
        # strip any formatting quotes and turn comma separated string into a list of fields
        gff_extra = parsed_args_be[0].gff_extra[0].strip("'\"").split(',')
    else:
        gff_extra = []

    # process the gff file to get required fields
    print("extra gff fields: " + str(gff_extra))
    gff_file = parsed_args_be[0].gff[0]
    gff_feat_type = ['CDS', 'tRNA', 'rRNA']

    min_depth_cutoff = parsed_args_be[0].min_depth[0]
    fraction_mismatch = parsed_args_be[0].mismatch[0]
    sam_file = parsed_args_be[0].bam[0]
    condition_label = parsed_args_be[0].label[0]
    print("\ncond label " + condition_label + "\n")

    # process pimms sam/bam  file and produce coordinate / bed files

    sam_dir, db_rdir, info_rdir = make_results_dirs_in_sam_dir(sam_file, condition_label)

    # process the gff file to get required fields
    gff_columns_addback, attr_to_columns = process_gff(gff_file, gff_feat_type, gff_extra, info_rdir)

    gff_columns_addback_pseudo, attr_to_columns_pseudo = process_gff(gff_file, ['pseudogene'], [], info_rdir)

    seqid_consistancy_check(gff_columns_addback, sam_file)
    process_sam(sam_file, min_depth_cutoff, fraction_mismatch)

    sam_stem = modify_sam_stem(sam_file, min_depth_cutoff, fraction_mismatch)

    # allocate insertions to features and create results merged with GFF
    # possibly poor coding to merge with gff here
    pimms_result_table_full = coordinates_to_features(sam_stem, attr_to_columns, gff_columns_addback, condition_label,
                                                      min_depth_cutoff, gff_extra, db_rdir)
    # if parsed_args[0].nano:
    #    print("--noreps forced for nanopore data\n")
    if not parsed_args_be[0].noreps:
        print("processing read groups as replicates for illumina data\n")
        mp_reps_feature_counts = coordinates_to_features_reps(sam_stem, attr_to_columns, condition_label)
        if not mp_reps_feature_counts.empty:
            merged_with_reps = pimms_result_table_full.merge(mp_reps_feature_counts, on=["seq_id", "start", "end"],
                                                             how='outer')
            # how='inner')
            pimms_result_table_full = merged_with_reps.fillna(0)
    else:
        print("not processing read groups as replicates\n")

    if not gff_columns_addback_pseudo.empty:
        tag_psueudogenes = gff_columns_addback_pseudo['locus_tag']
        pimms_result_table_full.loc[pimms_result_table_full.locus_tag.isin(tag_psueudogenes), "type"] = \
            pimms_result_table_full['type'] + '_pseudo'

    # print(parsed_args_be[0].out_fmt[0] + "out_fmt\n")
    # write results as text/excel
    if parsed_args_be[0].out_fmt[0] == 'tsv':
        pimms_result_table_full.to_csv(sam_stem + "_countinfo.tsv", index=False, sep='\t')
    elif parsed_args_be[0].out_fmt[0] == 'csv':
        pimms_result_table_full.to_csv(sam_stem + "_countinfo.csv", index=False, sep=',')
    else:
        writer = pd.ExcelWriter(sam_stem + '_countinfo.xlsx', engine='xlsxwriter')
        # Convert the dataframe to an XlsxWriter Excel object.
        pimms_result_table_full.to_excel(writer, sheet_name='PIMMS2_result', index=False)
        # Close the Pandas Excel writer and output the Excel file.
        writer.save()

    os.rename(sam_stem + "_countinfo." + parsed_args_be[0].out_fmt[0], os.path.join(db_rdir, sam_stem + "_countinfo." + parsed_args_be[0].out_fmt[0]))
    os.rename(sam_stem + "_insert_coords.txt", os.path.join(info_rdir, sam_stem + "_insert_coords.txt"))
    os.rename(sam_stem + ".bed", os.path.join(info_rdir, sam_stem + ".bed"))


# end bam_extract_func


def table_merge_func(parsed_args_tm):
    print(pimms_mssg + parsed_args_tm[0].command + pimms_mssg2)
    if parsed_args_tm[0].xlsx:
        print("Join: ", parsed_args_tm[0].xlsx[0], "\t", parsed_args_tm[0].xlsx[1], "\n")
        # requires  dependancy installed
        result_df1 = pd.read_excel(parsed_args_tm[0].xlsx[0], engine="openpyxl")
        result_df2 = pd.read_excel(parsed_args_tm[0].xlsx[1], engine="openpyxl")
        results_merged = pd.DataFrame.merge(result_df1, result_df2)
        writer = pd.ExcelWriter('merged_result.xlsx', engine='xlsxwriter')
        # Convert the dataframe to an XlsxWriter Excel object.
        results_merged.to_excel(writer, sheet_name='PIMMS2_merged_result', index=False)
        # Close the Pandas Excel writer and output the Excel file.
        writer.save()
    elif parsed_args_tm[0].csv:
        print("Join: ", parsed_args_tm[0].csv[0], "\t", parsed_args_tm[0].csv[1], "\n")
        result_df1 = pd.read_csv(parsed_args_tm[0].csv[0]).replace('"', '', regex=True)
        result_df2 = pd.read_csv(parsed_args_tm[0].csv[1]).replace('"', '', regex=True)
        results_merged = pd.DataFrame.merge(result_df1, result_df2)
        results_merged.to_csv('merged_result.csv', index=False)
    elif parsed_args_tm[0].tsv:
        print("Join: ", parsed_args_tm[0].tsv[0], "\t", parsed_args_tm[0].tsv[1], "\n")
        result_df1 = pd.read_csv(parsed_args_tm[0].tsv[0], sep="\t")
        result_df2 = pd.read_csv(parsed_args_tm[0].tsv[1], sep="\t")
        results_merged = pd.DataFrame.merge(result_df1, result_df2)
        results_merged.to_csv('merged_result.txt', index=False, sep="\t")
    else:
        print("\nUnable to merge results tables\n")


parsed_args = parse_arguments()  # parse command line arguments


# FIND_FLANK ###
def find_flank_func(parsed_args_ff):
    # if parsed_args[0].command == 'find_flank':
    print(pimms_mssg + parsed_args_ff[0].command + pimms_mssg2)
    # print((vars(parsed_args)))
    # sys.exit(1)
    # config_file = parsed_args.config_file[0]

    # construct config parser
    # p2config = configparser.ConfigParser()
    mapper = parsed_args_ff[0].mapper
    label = parsed_args_ff[0].label[0]
    print("\nFF label " + label + "\n")

    if parsed_args_ff[0].out_dir[0]:
        out_dir_ff = parsed_args_ff[0].out_dir[0]
    # print('\ncreating result dir: ' + out_dir + '\n')
    else:
        out_dir_ff = 'pimms2_' + label + '_' + time.strftime("%d%m%y_%H%M%S")
    # print('\ncreating result dir: ' + out_dir + '\n')
    # createFolder(out_dir)

    if os.path.isdir(out_dir_ff):
        print('\nresult dir exists\n')
    else:
        print('\ncreating result dir: ' + out_dir_ff + '\n')
        create_folder(out_dir_ff)

    out_dir_logs = os.path.join(out_dir_ff, 'logs')
    create_folder(out_dir_logs)

    fwdrev_wc = parsed_args_ff[0].fwdrev[0].strip("'\"").split(',')

    # exit(0)

    # print(pimms_mls)
    # dir(parsed_args[0].cpus[0])
    print('ncpus=' + str(parsed_args_ff[0].cpus[0]))
    # sys.exit(1)
    ncpus = int(parsed_args_ff[0].cpus[0])

    nano = parsed_args_ff[0].nano

    # experimental decontaminate transposon/vector sequence
    # not currently effective try another implementation when time allows?
    # decontam_tranposon = False
    #  print(parsed_args_ff[0].sub[0])
    fuzzy_levenshtein = bool(parsed_args_ff[0].lev[0])

    # set up some variables:
    if nano:  # nano == True
        # parsed_args[0].noreps = True
        fuzzy_levenshtein = True
        l_dist = parsed_args_ff[0].lev[0]  # maximum Levenshtein Distance
        min_length = 50
        max_length = 200
        print('overriding with Nanopore appropriate settings: Levenshtein distance of ' + str(
            l_dist) + ' + sequence length min = ' + str(min_length) + ', max = ' + str(max_length))
    else:
        subs = parsed_args_ff[0].sub[0]
        l_dist = parsed_args_ff[0].lev[0]  # maximum Levenshtein Distance
        insrt = parsed_args_ff[0].insert[0]
        dels = parsed_args_ff[0].deletion[0]
        min_length = parsed_args_ff[0].min
        max_length = parsed_args_ff[0].max

    # set up some names
    if nano:
        seqtype = '_nano'
    else:
        seqtype = ''

    if fuzzy_levenshtein:
        fq_result_suffix = (seqtype + "_pimms2out_trim" + str(max_length) + "_lev" + str(l_dist) + ".fastq")
    elif insrt > 0 | dels > 0:
        fq_result_suffix = (
                seqtype + "_pimms2out_trim" + str(max_length) + "_sub" + str(subs) + "_ins" + str(insrt) + "_del" + str(dels) + ".fastq")
    else:
        fq_result_suffix = (seqtype + "_pimms2out_trim" + str(max_length) + "_sub" + str(subs) + ".fastq")

    sam_result_suffix = re.sub('.fastq', '.sam', fq_result_suffix)

    fastq_dir = os.path.join(parsed_args_ff[0].in_dir[0], '')

    flanking_fastq_result_list = []

    if nano:  # nano == True

        glob_wc = ["*q.gz", "*.fastq", "*.fasta", "*.fasta.gz"]
        glob_read_files = find_read_files_with_glob(fastq_dir, glob_wc)
        print(glob_read_files, "...\n")
        print("nanopore PIMMS filtering starting...\n")
        print(datetime.datetime.now())

        pi = multiprocessing.Pool(ncpus)
        # for fq in glob.glob(fastq_dir + "*[aq].gz"):
        for fq in glob_read_files:
            fq_processed = os.path.join(out_dir_ff, Path(Path(fq).stem).stem + fq_result_suffix)
            flanking_fastq_result_list = flanking_fastq_result_list + [fq_processed]
            pi.apply_async(pimms_fastq,
                           # args=(fq, fq_processed, nano)
                           args=(fq, fq_processed, out_dir_logs, nano)
                           )

        pi.close()
        pi.join()

        # pi = multiprocessing.Pool(ncpus)
        # for fq in glob.glob(fastq_dir + "*.fast[aq]"):
        #     fq_processed = os.path.join(out_dir_ff, Path(Path(fq).stem).stem + fq_result_suffix)
        #     flanking_fastq_result_list = flanking_fastq_result_list + [fq_processed]
        #     pi.apply_async(pimms_fastq,
        #                    args=(fq, fq_processed)
        #                    )
        #
        # pi.close()
        # pi.join()

        print("nanopore PIMMS filtering completed...\n")
        print(datetime.datetime.now())

    else:  # nano == False

        glob_wc = ["*q.gz", "*.fastq"]
        glob_read_files = find_read_files_with_glob(fastq_dir, glob_wc)

        print("PIMMS read filtering starting...\n")
        print(datetime.datetime.now())
        pi = multiprocessing.Pool(ncpus)
        # for fq in glob.glob(fastq_dir + "*.fastq"):
        for fq in glob_read_files:
            if not (fwdrev_wc[0] in fq or fwdrev_wc[1] in fq):
                print("ERROR(fastq): text substrings " + fwdrev_wc[0] + "/" + fwdrev_wc[
                    1] + " NOT FOUND in read filenanes (to identify illumina fwd/rev fastq files)")
                print("ERROR(fastq): Check the fastq file names and/or update the --fwdrev parameter")
                sys.exit(1)
            fq_processed = os.path.join(out_dir_ff, Path(Path(fq).stem).stem + fq_result_suffix)
            pi.apply_async(pimms_fastq,
                           args=(fq, fq_processed, out_dir_logs, nano)
                           )

        pi.close()
        pi.join()

        # pi = multiprocessing.Pool(ncpus)
        # for fq in glob.glob(fastq_dir + "*q.gz"):
        #     if not (fwdrev_wc[0] in fq or fwdrev_wc[1] in fq):
        #         print("ERROR(fastq): text substrings " + fwdrev_wc[0] + "/" + fwdrev_wc[
        #             1] + " NOT FOUND in read filenanes (to identify illumina fwd/rev fastq files)")
        #         print("ERROR(fastq): Check the fastq file names and/or update the --fwdrev parameter")
        #         sys.exit(1)
        #
        #     fq_processed = os.path.join(out_dir_ff, Path(Path(fq).stem).stem + fq_result_suffix)
        #     pi.apply_async(pimms_fastq,
        #                    args=(fq, fq_processed, out_dir_logs)
        #                    )
        #
        # pi.close()
        # pi.join()
        print("PIMMS read filtering completed...\n")
        print(datetime.datetime.now())

        # match fwdrev match substrings e.g: _R1_/_R2_ --fwdrev parameter
        fqp_results_fwd = sorted(glob.glob(os.path.join(out_dir_ff, "*" + fwdrev_wc[0] + "*" + fq_result_suffix)))
        fqp_results_rev = sorted(glob.glob(os.path.join(out_dir_ff, "*" + fwdrev_wc[1] + "*" + fq_result_suffix)))
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

            tempfqname = Path(fwd_fqp_result).name
            # fix so only substring in file name nit dir is updated
            mrg_fqp_result_filename = os.path.join(Path(fwd_fqp_result).parent,
                                                   re.sub(fwdrev_wc[0], '_RX_', tempfqname,
                                                          count=1))  # replace  fwd substring '_R1_'
            # mrg_fqp_result_filename = re.sub(fwdrev_wc[0], '_RX_', fwd_fqp_result, count=1)
            flanking_fastq_result_list = flanking_fastq_result_list + [mrg_fqp_result_filename]
            print(mrg_fqp_result_filename)
            print(str(len(a)))
            print(str(len(b)))
            print(str(len(c)))
            # pysam bug doesn't parse files gzipped in chunks so gzipping removed here
            # with pysam.FastxFile(fwd_fqp_result) as fin, gzip.open(mrg_fqp_result_filename, mode='wt') as fout:
            with pysam.FastxFile(fwd_fqp_result, persist=False) as fin, open(mrg_fqp_result_filename,
                                                                             mode='wt') as fout:
                for entry in fin:
                    fout.write((str(entry) + '\n'))

            # with pysam.FastxFile(rev_fqp_result) as fin, gzip.open(mrg_fqp_result_filename, mode='at') as fout:
            with pysam.FastxFile(rev_fqp_result, persist=False) as fin, open(mrg_fqp_result_filename,
                                                                             mode='at') as fout:
                for entry in fin:
                    if entry.name in c:
                        fout.write((str(entry) + '\n'))
                        # fout.write(str(entry) + '\n')

        # remove intermediate fastq files
        # if parsed_args[0].rmfiles:
        if not parsed_args_ff[0].keep:
            delete_file_list(fqp_results_fwd)
            delete_file_list(fqp_results_rev)

        print("illumina merge of fwd/reverse data  completed...\n")

    # tidy up
    print(flanking_fastq_result_list)
    # concat_result_fastq = concat_fastq(flanking_fastq_result_list, parsed_args[0].label[0], fq_result_suffix, out_dir)
    concat_result_fastq = concat_fastq_raw(flanking_fastq_result_list, label, fq_result_suffix, out_dir_ff)

    # merge logs from different parallel cpus
    merge_logs(out_dir_logs)

    # do mapping stuff
    bam_name_ff = ''
    if parsed_args_ff[0].nomap:
        print("Skipping mapping step...\n")
    else:
        if mapper == 'minimap2' or nano:
            sam_output_mm = os.path.splitext(parsed_args_ff[0].fasta[0])[0] + '_' + label + re.sub('.sam', '_mm2.sam',
                                                                                                   sam_result_suffix)
            sam_output_mm = os.path.join(out_dir_ff, sam_output_mm)
            run_minimap2(concat_result_fastq, sam_output_mm, parsed_args_ff[0].fasta[0])
            bam_name_ff = py_sam_to_bam(sam_output_mm)
        elif mapper == 'bwa':
            sam_output_bwa = os.path.splitext(parsed_args_ff[0].fasta[0])[0] + '_' + label + re.sub('.sam', '_bwa.sam',
                                                                                                    sam_result_suffix)
            sam_output_bwa = os.path.join(out_dir_ff, sam_output_bwa)
            run_bwa(concat_result_fastq, sam_output_bwa, parsed_args_ff[0].fasta[0], ncpus)
            bam_name_ff = py_sam_to_bam(sam_output_bwa)

    return out_dir_ff, bam_name_ff


# end find_flank_func ###

# FIND_FLANK ###

if parsed_args[0].command == 'find_flank':
    out_dir, bam_name = find_flank_func(parsed_args)

# BAM_EXTRACT ###

elif parsed_args[0].command == 'bam_extract':
    bam_extract_func(parsed_args)


# TABLE_MERGE ###

elif parsed_args[0].command == 'table_merge':
    table_merge_func(parsed_args)


# FULL_PROCESS ###

elif parsed_args[0].command == 'full_process':
    out_dir, bam_name = find_flank_func(parsed_args)
    parsed_args[0].out_dir[0] = out_dir
    parsed_args[0].bam[0] = bam_name
    bam_extract_func(parsed_args)
