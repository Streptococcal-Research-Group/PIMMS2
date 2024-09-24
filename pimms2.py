import os, sys, subprocess, shutil
import datetime, time, re
from operator import itemgetter
import fileinput
import glob
import gzip
import multiprocessing
import random  # for log file names
import urllib as ul  # for removing url style encoding from gff text notes
from pathlib import Path
import configargparse
import pandas as pd
import pandasql as ps
import pysam  # sequence format specific module fastq/bam/sam...
import gffpandas.gffpandas as gffpd  # annotation format specific module gff3
from fuzzysearch import find_near_matches

pimms_mssg = """

\t===========================================================================================================
\tPragmatic Insertional Mutation Mapping system (PIMMS) mapping pipeline v2
\t===========================================================================================================
\t       o         o
\t       o         o
\t      //        //
\t     //        //
\t   |_||_|    |_||_|   @@@@@  @@@@@@  @@     @@  @@     @@   @@@@@@     @@@@@@
\t   |@||@|    |@||@|   @@  @@   @@    @@@@ @@@@  @@@@ @@@@  @@         @@    @@
\t   |@||@|    |@||@|   @@@@@    @@    @@ @@@ @@  @@ @@@ @@   @@@            @@
\t   |@||@|    |@||@|   @@       @@    @@  @  @@  @@  @  @@     @@@         @@
\t   |@@@@|    |@@@@|   @@       @@    @@     @@  @@     @@       @@      @@
\t   |@@@@|    |@@@@|   @@     @@@@@@  @@     @@  @@     @@  @@@@@@@    @@@@@@@@
\t===========================================================================================================
\tPIMMS2 mode: {0}
\t===========================================================================================================

"""


class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end


def create_folder(directory):
    try:
        os.makedirs(directory, exist_ok=True)
    except OSError:
        print('Error: Creating directory. ' + directory)


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


FIND_FLANK_MODE   = 'find_flank'
EXTRACT_BAM_MODE  = 'bam_extract'
FULL_PROCESS_MODE = 'full_process'
TABLE_MERGE       = 'table_merge'

ap = configargparse.ArgumentParser(
    prog="pimms2",
    add_config_file_help=False,
    config_file_parser_class=configargparse.DefaultConfigFileParser,
    epilog="\n\n*** N.B. This is a development version ***\n \n ",
    description='''description here''' )
    
ap.add_argument('-v', '--version', action='version', version='%(prog)s 2.1 demo')

modes = ap.add_subparsers(parser_class=configargparse.ArgParser, dest='command')

findflank = modes.add_parser(FIND_FLANK_MODE, 
                             add_config_file_help = False,
                             help                 = "Mode: find read regions flanking the IS sequence by mapping them to the target genome",
                             description          = "Args that start with '--' (eg. --fasta) can also be set in a config file (specified via -c)")

samcoords = modes.add_parser(EXTRACT_BAM_MODE, 
                             add_config_file_help = False,
                             help                 = "Mode: extract insertion site coordinates from bam file",
                             description          = "Args that start with '--' (eg. --fasta) can also be set in a config file (specified via -c)")

tablemerge = modes.add_parser(TABLE_MERGE, 
                              add_config_file_help = False,
                              help                 = 'Mode: merge two compatible PIMMS results tables '
                              '(N.B: this step does a simple table join and does not check the data)').add_mutually_exclusive_group()

fullprocess = modes.add_parser(FULL_PROCESS_MODE, 
                               add_config_file_help = False,
                               help                 = f"Mode: {FIND_FLANK_MODE} + {EXTRACT_BAM_MODE}",
                               description          = "Args that start with '--' (eg. --fasta) can also be set in a config file (specified via -c)")

# GENERAL ARGS

for subparser in [ findflank
                  ,samcoords
                  ,fullprocess]:

    subparser.add_argument("-c", "--config", 
                           required       = False, 
                           is_config_file = True,  # dest='config_file',
                           metavar        = 'pimms2.config',
                           help           = "use parameters from config file")

    subparser.add_argument("--nano", 
                           required = False, 
                           action   = 'store_true', 
                           default  = False,
                           help     = "override with settings more suitable for nanopore")

    subparser.add_argument("--label", 
                           required = False, 
                           metavar  = 'condition_name', 
                           default  = '',
                           help     = "identifying text tag to add to results file")


# FIND_FLANK ARGS
# to fix: nargs='?' deal with mistaken use of nargs=1 which give a single element list

for subparser in [ findflank
                  ,fullprocess ]:

    subparser.add_argument("--fasta", 
                           required = False, 
                           metavar  = 'ref_genome.fasta', 
                           type     = extant_file,
                           help     = "fasta file for reference genome ")

    subparser.add_argument("--qual_char", 
                           required = False, 
                           nargs    = '?', 
                           type     = str, 
                           default  = '0', 
                           choices  = [chr(x + 33) for x in list(range(12, 31))],
                           help     = "substitute a quality score ascii character when fasta read files used (nanopore only) (phred +33: ascii +:?) ['0']")

    subparser.add_argument("--nomap", 
                           required = False, 
                           action   = 'store_true', 
                           default  = False,
                           help     = "do not run mapping step")

    subparser.add_argument("--mapper", 
                           required = False, 
                           nargs    = '?', 
                           type     = str, 
                           default  = 'bwa', 
                           choices  = ['minimap2', 'bwa'],
                           help     = "select mapping software from available options")

    subparser.add_argument("--single", 
                           required = False, 
                           action   = 'store_true', 
                           default  = False,
                           help     = "only single direction Illumina data provided")

    subparser.add_argument("--keep", 
                           required = False, 
                           action   = 'store_true', 
                           default  = False,
                           help     = "keep intermediate fastq files etc for diagnostic purposes")

    subparser.add_argument("--lev", 
                           required = False, 
                           type     = int, 
                           default  = 0,
                           help     = "use Levenshtein distance (combined insert|del|sub score) [0]")

    subparser.add_argument("--sub", 
                           required = False, 
                           type     = int, 
                           default  = 0, # defaulted to 0 for consistency
                           help     = "number of permitted base substitutions in motif match [0]")

    subparser.add_argument("--insert", 
                           required = False, 
                           type     = int, 
                           default  = 0,
                           help     = "number of permitted base insertions in motif match [0]")

    subparser.add_argument("--del", 
                           required = False, 
                           type     = int, 
                           default  = 0, 
                           dest     = 'deletion',
                           help     = "number of permitted base deletions in motif match [0]")

    subparser.add_argument("--in_dir", 
                           required = True, 
                           dest     = 'in_dir', 
                           type     = extant_file,
                           help     = "directory containing input fastq files (assumed to match '*q.gz' or '*.fastq')")

    subparser.add_argument("--fwdrev", 
                           required = False, 
                           type     = str, 
                           default  = '_R1_,_R2_',
                           help     = "text substring to uniquely identify illumina fwd/rev paired fastq files ['_R1_,_R2_']")

    subparser.add_argument("--out_dir", 
                           required = False, 
                           metavar  = 'out_dir', 
                           default  = None,
                           action   = 'store',
                           help     = "directory to contain result files ['pimms2_`label`_`dmy`_`HMS`']")

    subparser.add_argument("--cpus", 
                           required = False, 
                           type     = int,
                           default  = int(os.cpu_count() / 2),
                           help     = "number of processors to use [(os.cpu_count() / 2)] ")

    subparser.add_argument("--max", 
                           required = False, 
                           type     = int, 
                           default  = 60,
                           help     = "clip results to this length [60]")

    subparser.add_argument("--min", 
                           required = False, 
                           type     = int, 
                           default  = 25,
                           help     = "minimum read length [25]")

    subparser.add_argument("--motif1", 
                           required = False, 
                           nargs    = 1, 
                           type     = str, 
                           default  = ['TCAGAAAACTTTGCAACAGAACC'],
                           # revcomp: GGTTCTGTTGCAAAGTTTTCTGA
                           help     = "IS end reference motif1 [TCAGAAAACTTTGCAACAGAACC](pGh9)")

    subparser.add_argument("--motif2", 
                           required = False, 
                           nargs    = 1, 
                           type     = str, 
                           default  = ['GGTTCTGTTGCAAAGTTTAAAAA'],
                           # revcomp: TTTTTAAACTTTGCAACAGAACC
                           help     = "IS end reference motif2 [GGTTCTGTTGCAAAGTTTAAAAA](pGh9)")


# BAM EXTRACT ARGS

for subparser in [ samcoords
                  ,fullprocess ]:

    subparser.add_argument("--bam", 
                           # exceptions for this option if bam extract only or full process
                           required = True if subparser == samcoords else False, 
                           metavar  = 'pimms.bam/sam', 
                           type     = extant_file,
                           help     = "bam/sam file of mapped IS flanking sequences " if subparser == samcoords else configargparse.SUPPRESS )

    ##@@ remove confusing mismatch parameter
    # subparser.add_argument("--mismatch",
    #                        required = False,
    #                        type     = float,
    #                        metavar  = 'float',
    #                        default  = None,
    #                        choices  = [round(x * 0.01, 2) for x in range(0, 21)],
    #                        help     = "fraction of permitted mismatches in mapped read ( 0 <= mismatch < 0.2) [no filter]")

    subparser.add_argument("--target_dup_offset",
                           required = False,
                           type     = int,
                           default  = 4,
                           metavar  = 'int',
                           help     = "target duplication offset (half duplication length)  int [4]")

    subparser.add_argument("--min_depth", 
                           required = False, 
                           type     = int, 
                           default  = 2, 
                           metavar  = 'int',
                           help     = "minimum read depth at insertion site >= int [2]")

    subparser.add_argument("--noreps", 
                           required = False, 
                           action   = 'store_true', 
                           default  = False,
                           help     = "do not separate illumina read groups as replicate insertion count columns")

    subparser.add_argument("--gff", 
                           required = True, 
                           type     = extant_file, 
                           metavar  = 'genome.gff',
                           help     = "GFF3 formatted file to use\n(note fasta sequence present in the file must be deleted before use)")

    subparser.add_argument("--gff_extra", 
                           required = False, 
                           type     = str, 
                           default  = [], 
                           metavar  = "'x,y,z'",
                           help     = "comma separated list of extra fields to include from the GFF3 annotation\ne.g. 'ID,translation,note' ")

    subparser.add_argument("--gff_force", 
                           required = False, 
                           action   = 'store_true', 
                           default  = False,
                           help     = "override GFF/BAM seq id discrepancies e.g. use when the gff has a plasmid not present in the reference sequence or vice-versa")

    subparser.add_argument("--out_fmt", 
                           required = False, 
                           type     = str, 
                           default  = 'xlsx',
                           choices  = ['xlsx', 'tsv', 'csv'],
                           help     = "set results table file format tab/comma separated or Excel (tsv|csv|xlsx) [xlsx]")


# TABLE_MERGE ARGS

tablemerge.add_argument("--xlsx", 
                        required = False, 
                        nargs    = 2, 
                        type     = extant_file,
                        help     = "2x .xlsx Excel files")

tablemerge.add_argument("--csv", 
                        required = False, 
                        nargs    = 2, 
                        type     = extant_file,
                        help     = "2x .csv comma separated text/table files")

tablemerge.add_argument("--tsv", 
                        required = False, 
                        nargs    = 2, 
                        type     = extant_file,
                        help     = '2x .tsv tab (\\t) separated text/table files')


parsed_args = ap.parse_known_args()

# check mode specified & at least one relevant flag/argument
if len(sys.argv) <= 2:
    # show usage message & exit
    ap.print_usage()
    sys.exit(1)

# show command line arguments specified
else:
    print("-----------------")
    ap.print_values()
    print("-----------------")


MODE_SELECTED = parsed_args[0].command

# show start message
print( pimms_mssg.format(MODE_SELECTED) )










######################
# FIND_FLANK COMMAND #
######################

def pimms_fastq(fqin_filename, fqout_filename ):

    # create complement dna translation matrix
    trans = str.maketrans('ATGCN', 'TACGN')

    # create motif type counts
    count = 0
    countqmulti = 0
    countqqrc = 0
    countq1q2 = 0
    wrongq2q1 = 0
    countq2rcq1rc = 0
    wrongq1rcq2rc = 0
    countq1 = 0
    countq1rc = 0
    countq2 = 0
    countq2rc = 0

    count_hits = { 'hit_q1_q2':     0, 'hit_but_short_q1_q2':     0,
                   'hit_q2rc_q1rc': 0, 'hit_but_short_q2rc_q1rc': 0,

                   'hit_q1_only':   0, 'hit_but_short_q1_only':   0,
                   'hit_q1rc_only': 0, 'hit_but_short_q1rc_only': 0,
                   
                   'hit_q2_only':   0, 'hit_but_short_q2_only':   0,
                   'hit_q2rc_only': 0, 'hit_but_short_q2rc_only': 0 }

    with pysam.FastxFile(fqin_filename, persist=False) as fin, open(fqout_filename, mode='wt') as fout:

        for entry in fin:

            # count fastq processed
            count += 1

            # specify transposon motif query; forward
            qry1 = parsed_args[0].motif1[0].strip("'\"")
            qry2 = parsed_args[0].motif2[0].strip("'\"")

            # specify transposon motif query; reverse complement
            qry1rc = qry1.translate(trans)[::-1]
            qry2rc = qry2.translate(trans)[::-1]
            
            # find transposon motif matches for each query
            matchseqs = [ find_near_matches(qry, 
                                            entry.sequence,                         # function defaults are None
                                            max_substitutions = SUBSTITUTIONS if not FUZZY_LEVENSHTEIN else None,
                                            max_deletions     = DELETIONS     if not FUZZY_LEVENSHTEIN else None,
                                            max_insertions    = INSERTIONS    if not FUZZY_LEVENSHTEIN else None,
                                            max_l_dist        = L_DIST        if     FUZZY_LEVENSHTEIN else None )
                          for qry in [ qry1
                                      ,qry2
                                      ,qry1rc
                                      ,qry2rc ] ]
            
            # unpack query matches
            matchesq1, matchesq2, matchesq1rc, matchesq2rc = matchseqs
            
            # check matches found
            if not any(matchseqs):
                continue

            # check not multiple matches for same motif
            if max( len(matches) for matches in matchseqs )  > 1:
                countqmulti += 1
                continue
            
            # specify further multiple motif match criteria
            multiple_qqrc_matchseqs = [
                # multiple matches to same motif ( query & reverse complement )
                ((len(matchesq1) == 1) and (len(matchesq1rc) == 1)),
                ((len(matchesq2) == 1) and (len(matchesq2rc) == 1)),

                # multiple matches to two incompatible motifs
                ((len(matchesq1) == 1) and (len(matchesq2rc) == 1)),
                ((len(matchesq2) == 1) and (len(matchesq1rc) == 1)) ]

            # check not further multiple motif matches
            if any(multiple_qqrc_matchseqs):
                countqqrc += 1
                continue



            # process motif match pairs to extract target sequences
            if (len(matchesq1) == 1) and (len(matchesq2) == 1):
                countq1q2    += 1
                capture_start = matchesq1[0].end
                capture_stop  = matchesq2[0].start
                hit_key       = 'hit_q1_q2'
                short_key     = 'hit_but_short_q1_q2' 

                if matchesq2[0].start <= matchesq1[0].end:
                    wrongq2q1 += 1
                    continue

                else:
                    report_fq     = True
                    report_rc     = False

            elif (len(matchesq1rc) == 1) and (len(matchesq2rc) == 1):
                countq2rcq1rc += 1
                capture_start  = matchesq2rc[0].end
                capture_stop   = matchesq1rc[0].start
                hit_key        = 'hit_q2rc_q1rc'
                short_key      = 'hit_but_short_q2rc_q1rc'

                if matchesq1rc[0].start <= matchesq2rc[0].end:
                    wrongq1rcq2rc += 1
                    continue

                else:
                    report_fq      = True
                    report_rc      = False


            # process single motif matches to extract target sequences
            elif len(matchesq1) == 1:
                countq1      += 1
                capture_start = matchesq1[0].end
                capture_stop  = None  # until end
                hit_key       = 'hit_q1_only'
                short_key     = 'hit_but_short_q1_only'
                report_fq     = True
                report_rc     = False


            elif len(matchesq2rc) == 1:
                countq2rc    += 1
                capture_start = matchesq2rc[0].end
                capture_stop  = None  # until end
                hit_key       = 'hit_q2rc_only'
                short_key     = 'hit_but_short_q2rc_only'
                report_fq     = True
                report_rc     = False


            elif len(matchesq1rc) == 1:
                countq1rc    += 1
                capture_start = 0
                capture_stop  = matchesq1rc[0].start
                hit_key       = 'hit_q1rc_only'
                short_key     = 'hit_but_short_q1rc_only'
                report_fq     = True
                report_rc     = True

            elif len(matchesq2) == 1:
                countq2      += 1
                capture_start = 0
                capture_stop  = matchesq2[0].start
                hit_key       = 'hit_q2_only'
                short_key     = 'hit_but_short_q2_only'
                report_fq     = True
                report_rc     = True

            # to catch any match combination not specified above
            else:
                report_fq     = False


            # report motif matches
            if report_fq:

                # extract motif match sequence info
                captured_seqstring  = str(entry.sequence)[capture_start:capture_stop]
                captured_qualstring = str(entry.quality )[capture_start:capture_stop]

                # adjust sequence quality scores for nanopore reads as required
                if NANOPORE_MODE:

                    if len(captured_qualstring) < 5:

                        captured_qualstring = QUAL_CHAR * len(captured_seqstring)

                # record short hit
                if not len(captured_seqstring) >= MIN_LENGTH:

                    count_hits[short_key] += 1

                # record hit
                else:

                    count_hits[hit_key] += 1

                    # trim match sequence; as forward
                    if not report_rc:

                        trimmed_seqstring  = captured_seqstring[0:MAX_LENGTH]
                        trimmed_qualstring = captured_qualstring[0:MAX_LENGTH]

                    # trim match sequence; as reverse complement
                    else:

                        trimmed_seqstring = captured_seqstring[-MAX_LENGTH:].translate(trans)[::-1]
                        trimmed_qualstring = captured_qualstring[-MAX_LENGTH:][::-1]

                    # ouput fastq
                    fout.write(f'@{str(entry.name)} CO:Z:{str(entry.comment)}\n')
                    fout.write(f'{trimmed_seqstring}\n')
                    fout.write('+\n')
                    fout.write(f'{trimmed_qualstring}\n')

    # specify log file
    log_file = f'{LOG_DIR}/log_{os.getppid()}_{multiprocessing.current_process().pid}_{random.randint(1000, 9999)}.txt'

    # open log & record summary
    with open( log_file, mode='w' ) as log:

        try:

            # extract relevant hit counts
            both_motif_total_counts = itemgetter( 'hit_q1_q2',
                                                  'hit_q2rc_q1rc',
                                                  'hit_but_short_q1_q2',
                                                  'hit_but_short_q2rc_q1rc' )(count_hits)

            both_motif_passed_counts = itemgetter( 'hit_q1_q2',
                                                   'hit_q2rc_q1rc' )(count_hits)

            sing_motif_passed_counts = itemgetter( 'hit_q1_only',
                                                   'hit_q2_only', 
                                                   'hit_q1rc_only', 
                                                   'hit_q2rc_only' )(count_hits)

            total_motif_passed_counts = both_motif_passed_counts + sing_motif_passed_counts

            # specify relevant headers / counts
            log_info = [ ( 'fastq',                              os.path.basename(fqin_filename)  ),
                         ( 'read count',                         count                            ),
                         ( 'multimatched',                       countqmulti                      ),
                         ( 'mismatched',                         countqqrc                        ),
                         ( 'fwd|rev matched',                    sum( both_motif_total_counts   ) ),
                         ( f'fwd|rev matched >= {MIN_LENGTH}bp', sum( both_motif_passed_counts  ) ),
                         ( f'single matched >= {MIN_LENGTH}bp',  sum( sing_motif_passed_counts  ) ),
                         ( 'total matched',                      sum( total_motif_passed_counts ) ) ]
        
            # extract headers / counts
            log_info = zip( *log_info )

            for info in log_info:

                print( *info, sep='\t', file=log )

        except Exception as e:

            log.write(f'ERROR: {e} problem with motif matching to {fqin_filename}!\n')

        # log closes automatically at 'with' 'as' end


# define general mode arguments
if MODE_SELECTED in [FIND_FLANK_MODE, EXTRACT_BAM_MODE, FULL_PROCESS_MODE]:

    NANOPORE_MODE    = parsed_args[0].nano
    RESULTS_LABEL    = parsed_args[0].label



# run find flank functions
if MODE_SELECTED in [FIND_FLANK_MODE, FULL_PROCESS_MODE]:

    # define specific find flank mode arguments
    NCPUS             = parsed_args[0].cpus
    
    IN_DIR            = parsed_args[0].in_dir
    OUT_DIR           = parsed_args[0].out_dir

    L_DIST            = parsed_args[0].lev
    MIN_LENGTH        = parsed_args[0].min
    MAX_LENGTH        = parsed_args[0].max

    SINGLE_END_MODE   = parsed_args[0].single

    FWDREV_IDS        = parsed_args[0].fwdrev
    
    SKIP_MAPPING     = parsed_args[0].nomap
    REFERENCE_PATH   = parsed_args[0].fasta
    MAPPING_SOFTWARE = parsed_args[0].mapper

    KEEP_ALL_FILES   = parsed_args[0].keep

    ##### do we also need to overide L_DIST to 1 for NANOPORE_MODE? currently set as default (0)
    FUZZY_LEVENSHTEIN = True if NANOPORE_MODE else bool(L_DIST)
    #####

    QUAL_CHAR         = parsed_args[0].qual_char

    SUBSTITUTIONS     = parsed_args[0].sub
    INSERTIONS        = parsed_args[0].insert
    DELETIONS         = parsed_args[0].deletion


    # append '/' to input directory as required
    FASTQ_DIR = os.path.join(IN_DIR, '')

    # specify output / log directory
    if not OUT_DIR:
        
        OUT_DIR = f'./pimms2'
        
        if RESULTS_LABEL:
            
            OUT_DIR += f'_{RESULTS_LABEL}'
        
        OUT_DIR += f'_{time.strftime("%d%m%y_%H%M%S")}'

    LOG_DIR = os.path.join(OUT_DIR, 'logs')


    # specify empty file suffix
    result_suffix = ''

    # append nanopore identifer
    if NANOPORE_MODE: 
        
        result_suffix += '_nano'

    # append max_length info
    result_suffix += f'_pimms2out_trim{MAX_LENGTH}'

    # append fuzzy match threshold info
    if FUZZY_LEVENSHTEIN:

        result_suffix += f'_lev{L_DIST}'

    else:
        
        # include substitution threshold info
        result_suffix += f'_sub{SUBSTITUTIONS}'

        # include insert/deletion threshold info
        if INSERTIONS > 0 | DELETIONS > 0:
    
            result_suffix += f'_ins{INSERTIONS}_del{DELETIONS}'
    
    # append relevant file extensions
    result_suffixes = [ f'{result_suffix}.{extension}'
                          for extension in ( 'fastq', 'sam', 'bam' ) ]

    # extract final fastq/sam file suffixes
    fq_result_suffix, sam_result_suffix, bam_result_suffix = result_suffixes




# PRE-RUN CHECKS

    # MAPPING CHECKS

    if not SKIP_MAPPING:
        
        # inform user that selected mapping software not installed
        if not shutil.which(MAPPING_SOFTWARE):

            print(f'ERROR: {MAPPING_SOFTWARE.upper()} not found in $PATH.\n')
            print(f'EXITING: install {MAPPING_SOFTWARE.upper()} or select another mapping software.\n')
            sys.exit(0)

        # inform user that reference genome not specified
        if not REFERENCE_PATH:
            print(f'ERROR: Reference genome not specified.\n')
            print(f'EXITING: Specify "--fasta /path/to/reference.fasta" or select "--nomap".\n')
            sys.exit(0)

        # inform user that results label not specified
        if not RESULTS_LABEL:

            print(f'ERROR: Results label not specified.\n')
            print(f'EXITING: Specify "--label LABEL_NAME" or select "--nomap".\n')
            sys.exit(0)    



    # INPUT FILES CHECK

    # split comma seperate paired read ids
    FWDREV_IDS = FWDREV_IDS.strip("'\"").split(',')

    # check enough paired read ids
    if len(FWDREV_IDS) != 2:

        print(f'ERROR: Inadequate paired reads IDs provided\n')
        print(f'EXITING: Specify two IDs delimited by a comma i.e. "--fwdrev_wc _ID1_,_ID2_".\n')
        sys.exit(0)
    
    # extract individual foward / reverse read ids
    else:

        FWD_ID, REV_ID = FWDREV_IDS


    # specify general search wild cards for fastqs
    glob_wc = [ "*q.gz", "*.fastq" ]
    
    # extend search wild cards to include additional nanopore specific fastq suffixes
    if NANOPORE_MODE:

        glob_wc.extend([ "*.fasta", "*.fasta.gz" ])

    # create store for found wildcard matches
    glob_read_files = [] 

    # cycle through search wild cards    
    for suffix_wc in glob_wc:

        # search fastq directory for wildcard matches
        wc_matches = glob.glob(f'{FASTQ_DIR}{suffix_wc}')

        # record any wildcard matches in store
        glob_read_files.extend(wc_matches)

    # check fastqs found
    if not glob_read_files:
    
        print(f'ERROR: No fastqs found in {FASTQ_DIR}.\n')
        print(f'EXITING: Ensure files have the following suffixes; {" ".join(glob_wc)}\n')
        sys.exit(0)

    # check fastq paired read id info present 
    if not NANOPORE_MODE and not SINGLE_END_MODE:                

        paired_ids_found = [ (FWD_ID in fq or REV_ID in fq)
                             for fq in glob_read_files ]

        if not all(paired_ids_found):

            print(f'ERROR: Not all fastq files contain the paired read IDs "{FWD_ID}" or "{REV_ID}".')
            print('EXITING: Rename fastq files and/or provide different IDs i.e. "--fwdrev_wc _ID1_,_ID2_".')
            sys.exit(0)




# DISPLAY RUN SETTINGS

    # show specified settings
    print(f'\n\t{MODE_SELECTED.upper()} SETTINGS:\n')
    print(f'\tNumber CPUs:          {NCPUS}')
    print(f'\tFastq directory:      {IN_DIR}')
    print(f'\tOutput directory:     {OUT_DIR}')
    print(f'\tLevenshtein distance: {L_DIST}')
    print(f'\tFuzzy matching:       {FUZZY_LEVENSHTEIN}')
    print(f'\tMinimum read length:  {MIN_LENGTH}')
    print(f'\tMaximum read length:  {MAX_LENGTH}')
    print(f'\tForward read ID:      {FWD_ID}')
    print(f'\tReverse read ID:      {REV_ID}')
    print(f'\tResults label:        {RESULTS_LABEL}')    
    

    READ_MODE = 'nanopore' if NANOPORE_MODE else 'illumina'
    
    # show read type settings
    print(f'\n\n\tREAD TYPE SETTINGS:\n')
    print(f'\tRead type:            {READ_MODE}')

    # show nanopore specific settings
    if NANOPORE_MODE: 

        print(f'\tQuality character:    {QUAL_CHAR}')

    # show illumina specific settings
    else:
    
        print(f'\tSubstitutions:        {SUBSTITUTIONS}')
        print(f'\tInsertions:           {INSERTIONS}')
        print(f'\tDeletions:            {DELETIONS}')

    # show mapping specific settings
    if not SKIP_MAPPING:

        print(f'\n\n\tMAPPING SETTINGS:\n')
        print(f'\tReference genome:     {REFERENCE_PATH}')
        print(f'\tMapping Software:     {MAPPING_SOFTWARE}')
        print(f'\tInstallation used:    $PATH/{MAPPING_SOFTWARE}')
    
    print('\n')



# PREPARE OUTPUTS DIRECTORIES

    # create output directory
    if not os.path.isdir(OUT_DIR):
        create_folder(OUT_DIR)

    # remove output logs from previous runs
    if os.path.isdir(LOG_DIR):  
        shutil.rmtree(LOG_DIR)

    # create output log directory
    if not os.path.isdir(LOG_DIR):
        create_folder(LOG_DIR)




# FILTER READS

    print(f'\n{len(glob_read_files)} input files found:\n')
    
    for fq in sorted(glob_read_files):

        fq = os.path.basename(fq)
        print(f'- {fq}')

    print(f'\n\nPIMMS read filtering started at {datetime.datetime.now()}...\n')

    pi = multiprocessing.Pool(NCPUS)

    for fq in glob_read_files:

        ##### as is likely lead to errors if only 1 suffix & contains (.)?; explicitly remove collective suffixes

        # extract fastq basename
        fq_basename = Path(fq).stem

        # extract fastq filename i.e. remove fastq suffixes
        fq_filename = Path(fq_basename).stem

        # specify fastq path
        fq_processed = os.path.join(OUT_DIR, f'{fq_filename}{fq_result_suffix}')

        print(f'Processing {os.path.basename(fq)} ==> {os.path.basename(fq_processed)}')

        # process fastq?
        pi.apply_async( pimms_fastq, args=(fq, fq_processed ) )

    pi.close()
    pi.join()

    print(f"\nPIMMS read filtering completed at {datetime.datetime.now()}...\n")



# PROCESS FILTERED OUTPUTS

    flanking_fastq_result_list = []

    # find filtered fastqs

    # find filtered files
    if NANOPORE_MODE or SINGLE_END_MODE:

        # specify filtered fastq search wild card
        fqs_filtered_wc = os.path.join(OUT_DIR, f'*{fq_result_suffix}')
    
        # search output directory for filtered fastq
        fqs_results = sorted( glob.glob(fqs_filtered_wc) )
        
        # store filtered fastq paths
        flanking_fastq_result_list.extend(fqs_results)


    # find filtered files & merge paired R1/R2 fastqs from same lanes
    else:

        print(f"\nMerging paired Illumina reads started at {datetime.datetime.now()}...\n")

        # specify filtered fastq search wild card
        fqp_filtered_wc = [ os.path.join(OUT_DIR, f'*{ID}*{fq_result_suffix}')
                            for ID in [ FWD_ID
                                        ,REV_ID ] ]

        fqp_results = [ sorted(glob.glob(fq_wc))
                        for fq_wc in fqp_filtered_wc ]
        
        fqp_results_fwd, fqp_results_rev = fqp_results

        ##### shouldn't assume 1 & 2 mathces found for all i.e. zipped correctly
        ##### if no reads pass filtering then filtered fastq won't exist
        for paired_results in zip(fqp_results_fwd, fqp_results_rev):
        #####
            fwd_fqp_result, rev_fqp_result = paired_results

            # specify survey fastq info
            survey_info = { 1: { 'fqout': fwd_fqp_result },

                            2: { 'fqout': rev_fqp_result } }

            # cycle through paired read members
            for pair_member in survey_info:

                # extract pair member sub dictionary
                member_info = survey_info[pair_member]

                # open pair member filtered fastq
                with pysam.FastxFile( member_info['fqout'], persist=False ) as fh:

                    # extract extracted names, remove duplicates & store
                    member_info['names'] = set( entry.name 
                                                for entry in fh )

            # unpack unique read names
            result_reads_unique = [ survey_info[pair_member]['names']
                                    for pair_member in sorted( survey_info ) ]

            # extract unique read name sets
            result1_reads_unique, result2_reads_unique = result_reads_unique
            
            # extract unique reverse reads i.e. no forward specific names & no forward/reverse shared names
            reads_difference = result2_reads_unique.difference(result1_reads_unique)

            # extract filtered fastq basenames
            fltrd_fq_basenames = [ os.path.basename(path) 
                                    for path in paired_results]

            fwd_fq_basename, rev_fq_basename = fltrd_fq_basenames

            # extract common fastq prefix
            fltrd_fq_prefix = os.path.commonprefix(fltrd_fq_basenames)
            
            # extract common suffix (invert, extract "prefix" then revert)
            fltrd_fq_suffix = os.path.commonprefix([ name[::-1] for name in fltrd_fq_basenames])[::-1]
            
            # specify merged fastq basename using common prefix/sufix but with distinct ID (X)
            merged_fq_basename = f'{fltrd_fq_prefix}X{fltrd_fq_suffix}'
            
            # specify merged fastq path
            merged_fqp_result_path = os.path.join( OUT_DIR,
                                                    merged_fq_basename )

            # store merged fastq results
            flanking_fastq_result_list.append( merged_fqp_result_path )

            # pysam bug doesn't parse files gzipped in chunks so gzipping removed here
            
            ##### as discussed am a bit unsure about logic here. 
            ##### also should shouldn't lanes be merged and not reads? all gets concatenated later anyhow
            print(f'Merging {fwd_fq_basename} + {rev_fq_basename} ==> {merged_fq_basename}')

            with open(merged_fqp_result_path, mode='wt') as fout:

                with pysam.FastxFile(fwd_fqp_result, persist=False) as fqfwd:

                    for entry in fqfwd:
                
                        # output any forward reads (if unique to forward fastq & shared between paired fastqs)
                        fout.write(f'{entry}\n')

                with pysam.FastxFile(rev_fqp_result, persist=False) as fqrev:

                    for entry in fqrev:

                        # output reverse reads only (if unique to reverse fastq & NOT shared between paired fastqs)
                        if entry.name in reads_difference:

                            fout.write(f'{entry}\n')
        
        print(f"\nMerging paired Illumina reads completed at {datetime.datetime.now()}...\n")


        # remove intermediate fastqs
        if not KEEP_ALL_FILES:

            print(f"\nRemoving intermediate filtered paired reads\n")

            for intermediates in fqp_results:

                delete_file_list(intermediates)



    # concatenate filtered fastqs

    # specify concatenated fastq path
    concat_result_fastq = os.path.join(OUT_DIR,
                                       f'{RESULTS_LABEL}_RX_concat{fq_result_suffix}.gz')

    ##### probably just use zcat; python module overcomplicating things?
    with gzip.open(concat_result_fastq, "wt", compresslevel=6) as big_file:
    #####

        with fileinput.input( files=flanking_fastq_result_list ) as inputs:
        
            for line in inputs:

                big_file.write(line)

    if not KEEP_ALL_FILES:

        print('\nRemoving intermediate fastq flanking reads files\n')
        
        delete_file_list(flanking_fastq_result_list)


    # merge logs from different parallel cpus

    # specify log file search wildcard
    log_wc = os.path.join( LOG_DIR,
                          "log_*txt" )

    # search for log files
    log_files = glob.glob( log_wc )

    # reads log files as dataframe
    df_from_each_log = ( pd.read_table(log) 
                         for log in log_files )

    # concatenate log files
    merged_logs_df = pd.concat( df_from_each_log, 
                                ignore_index=True )

    # sort logs by fastq filename
    merged_logs_df = merged_logs_df.sort_values( by=['fastq'] )
    
    # calculate metric totals
    log_sums = merged_logs_df.sum(numeric_only=True)
    
    log_sums['fastq'] = 'COMBINED'
    
    merged_logs_df = merged_logs_df.append(log_sums, ignore_index=True)
    
    # specify merged log path
    merged_log_path = os.path.join( OUT_DIR, 
                                    'result_summary.txt' )
    
    # write merged log
    merged_logs_df.to_csv( merged_log_path, 
                           sep='\t', 
                           index=False )
    
    # print filtered metrics to screen

    print('\nMOTIF MATCHES:\n')

    print(f'{merged_logs_df.to_string(index=False)}\n')




# MAP READS


    if SKIP_MAPPING: 
        
        print("Skipping mapping step...\n")

        bam_name = None

    else:

        # extract reference directory / filename
        ref_dir, ref_file = os.path.split( REFERENCE_PATH )

        # extract reference basename
        ref_basename, *_ = os.path.splitext( ref_file )

        software_info = { 'minimap2': 'mm2'
                         ,'bwa':      'bwa' }

        software_id = software_info[MAPPING_SOFTWARE]
        
        mapped_files = [ f'{ref_basename}_{software_id}_{RESULTS_LABEL}{suffix}'
                         for suffix in (sam_result_suffix, bam_result_suffix) ]

        mapped_paths = [ os.path.join( OUT_DIR, name )   
                         for name in mapped_files ]

        sam_path, bam_path = mapped_paths

        if MAPPING_SOFTWARE == 'minimap2' or NANOPORE_MODE:
            
            # specify mapping mode
            mm_mode = 'map-ont' if NANOPORE_MODE else 'sr'

            # specify minimap2 command       
            map_cmd = f'minimap2 -x {mm_mode} -a -y -o {sam_path} {REFERENCE_PATH} {concat_result_fastq} --secondary=no --sam-hit-only'


        elif MAPPING_SOFTWARE == 'bwa':

            idx_dir = f'{ref_basename}_index'

            # specify index directory reference path
            idx_ref_path = os.path.join( idx_dir,
                                         ref_file )

            idx_sa_path = os.path.join( idx_dir,
                                        f'{ref_file}.sa' )
            
            if os.path.exists(idx_sa_path):
                
                print('Using existing BWA index...')
            
            else:

                print('Creating BWA index...')

                create_folder(idx_dir)
                
                # copy reference fasta to index directory
                shutil.copyfile( REFERENCE_PATH, 
                                 idx_ref_path)

                # specify bwa index command
                idx_cmd = f'bwa index {idx_ref_path}'

                # genereate index
                print(f'RUNNING: {idx_cmd}')
                process = subprocess.Popen( idx_cmd.split(' '), 
                                            stdout=subprocess.PIPE, 
                                            stderr=subprocess.PIPE )
                stdout, stderr = process.communicate()
                print(stdout.decode('utf-8'))
                print(stderr.decode('utf-8'))

            # specify bwa command
            map_cmd = f'bwa mem -t {NCPUS} -C -o {sam_path} {idx_ref_path} {concat_result_fastq}'


        # perform mapping
        print(f'RUNNING: {map_cmd}')
        process = subprocess.Popen( map_cmd.split(' '), 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE )
        stdout, stderr = process.communicate()
        print(stdout.decode('utf-8'))
        print(stderr.decode('utf-8'))


        # CONVERT SAM -> BAM

        # add?? line to remove unmapped readsbased on: | samtools view -F 4 -o onlyMapped.bam ??
        # noinspection PyUnresolvedReferences
        pysam.sort( '-O', 'BAM', "-o", bam_path, sam_path )

        # noinspection PyUnresolvedReferences
        pysam.index( bam_path )

        print('\nMapping stats (flagstat):\n')

        # noinspection PyUnresolvedReferences
        for fsline in pysam.flagstat(bam_path).splitlines()[:5]:

            print(fsline)
        
        # remove intermediate sam
        if not KEEP_ALL_FILES:

            delete_file_list([sam_path])










#######################
# BAM_EXTRACT COMMAND #
#######################

# run bam extract functions
if MODE_SELECTED in [ EXTRACT_BAM_MODE, FULL_PROCESS_MODE ]:
        

    # define specific bam extract mode arguments

    MIN_DEPTH         = parsed_args[0].min_depth

    TARGET_DUP_OFFSET = parsed_args[0].target_dup_offset
    ##@@ FRACTION_MISMATCH = parsed_args[0].mism
    # BAM_PATH          = parsed_args[0].bam    if EXTRACT_BAM_MODE else bam_path
    BAM_PATH          = parsed_args[0].bam    if MODE_SELECTED == EXTRACT_BAM_MODE else bam_path

    # print(f'BAM_PATH: {BAM_PATH}')
    # print(f'bam_path: {bam_path}')
    # print(f'MODE_SELECTED: {MODE_SELECTED}')

    SKIP_REPS         = parsed_args[0].noreps if not NANOPORE_MODE else True

    GFF_PATH          = parsed_args[0].gff

    GFF_EXTRA         = parsed_args[0].gff_extra
    
    GFF_FORCE         = parsed_args[0].gff_force
    
    OUT_FORMAT        = parsed_args[0].out_fmt
    
    # convert fields to list
    if GFF_EXTRA:

        GFF_EXTRA = GFF_EXTRA.strip("'\"").split(',')  

    OUT_DIR = os.path.dirname( BAM_PATH )
    
    # N.B. same as outdir specified in find_flank mode

    GFF_name = os.path.basename( GFF_PATH )

    GFF_basename, *_ = os.path.splitext( GFF_name )


    # extract alignment file name
    BAM_basename = os.path.basename( BAM_PATH )

    # remove alignment file extension
    BAM_filename, *_ = os.path.splitext( BAM_basename )

    # append depth & mismatch info
    ##@@ BAM_filename += f'_md{MIN_DEPTH}_mm{FRACTION_MISMATCH or 0}'
    BAM_filename += f'_md{MIN_DEPTH}'



    # specify out subdirectory labels
    out_subdir_labels    = ( 'out_dashboard','out_info' )

    # specify results subdirectory labels
    timeinfo = time.strftime("%d%m%y_%H%M%S")
    result_subdir_labels = (f'{timeinfo}_results_dashboard', f'{timeinfo}_results_info' )

    # specify output subdirectory paths
    subdirs = [ [ os.path.join(OUT_DIR, f'{RESULTS_LABEL}_{label}' ) for label in labels ] 

                  for labels in [ out_subdir_labels, result_subdir_labels ] ]
    
    # unpack different subdirs
    (out_dash_subdir, *_), *_ = out_subdirs, results_subdirs = subdirs
    # N.B. first subdirectory used to assess whether to use out or results subdirectories

    # select relevant output subdirectories
    relelvant_subdirs = out_subdirs if not os.path.exists( out_dash_subdir ) else results_subdirs

    # unpack individual relevant subdirectories
    db_rdir, info_rdir = relelvant_subdirs
    
    # create relevant output subdirectories
    try:
        [ os.makedirs( subdir, exist_ok=True ) for subdir in relelvant_subdirs ]
    
    except OSError:
        print(f'Error while creating result dirs in {OUT_DIR}')



    # specify output file names
    outputs = [ os.path.join( subdirectory, 
                              f'{BAM_filename}{suffix}' )
                for subdirectory, suffix in ( ( info_rdir, '.bed'                    ),
                                              ( info_rdir, '_insert_coords.txt'      ),
                                              ( db_rdir,   f'_countinfo.{OUT_FORMAT}') ) ]
    
    # specify countinfo output file
    bed_file, coordinate_file, countinfo_file = outputs



    # process gff
    
    GFF_info = {

        1: { 'features':   ['CDS', 'tRNA', 'rRNA'],
             'extras':     [*GFF_EXTRA],
             'processed':  [],
             'pseudogene': False },

        2: { 'features':   ['pseudogene'],
             'extras':     [],
             'processed':  [],
             'pseudogene': True } }


    # specify standard gff annotation to report
    gff_standard_fields = [ 'seq_id'      # standard
                           ,'type'        # standard
                           ,'start'       # standard
                           ,'end'         # standard
                    
                           ,'feat_length' # added
                           ,'product'     # added

                           ,'ID'          # attribute
                           ,'locus_tag'   # attribute
                           ,'gene' ]      # attribute
                
    ##### ATTRIBUTES ARE INCONSISTENT ACCROSS GFFS! CANNOT GUARANTEE THEY EXIST


    for key in GFF_info:

        # specify gff identifier
        GFF_filtered_identifier = '_pseudogene' if GFF_info[key]['pseudogene'] else ''

        # specify filtered gff filename
        GFF_filtered_name = f'{GFF_basename}_pimms_features{GFF_filtered_identifier}.gff'

        # specify filtered
        GFF_filtered_path = os.path.join( info_rdir, GFF_filtered_name )

        # read in gff annotation
        GFF_annotation = gffpd.read_gff3(GFF_PATH)

        # filter gff annotation for relevant features
        GFF_annotation = GFF_annotation.filter_feature_of_type( GFF_info[key]['features'] )

        # write filtered gff to output        
        GFF_annotation.to_gff3( GFF_filtered_path )

        # convert filtered annotation to dataframe & split final attributes column into individual columns 
        GFF_df = GFF_annotation.attributes_to_columns()
        # N.B. indivudal attributes recorded only if relevant i.e. cell values if exists for entry otherwise 'None'

        # filtering removed all rows (no relevant features)
        if GFF_df.empty:

            # columns; seq_id, source, type, start, end, score, strand, phase, attributes
            GFF_info[key]['processed'].extend([GFF_df, GFF_df])
        
        else:

            # create new column for feature length; difference between feature start & feature end
            GFF_df = GFF_df.assign( feat_length=(GFF_df.end - GFF_df.start + 1) )

            # remove columns if all blank values
            # remove original attributes column; split columns for individual attributes retained
            GFF_df = GFF_df.dropna( axis=1
                                   ,how='all' ).drop(columns=['attributes'])

            # create product column if absent
            if not 'product' in GFF_df.columns:
            
                GFF_df['product'] = '-'
            
            else:
                # you've just removed column if na so why fill?
                GFF_df['product'] = GFF_df['product'].fillna('').astype(str).apply(ul.parse.unquote)  # added fix for None datatype


            # fix to skip requested extra gff annotation field if not present in GFF
            attributes_absent = [ attribute 
                                  for attribute in GFF_info[key]['extras'] 
                                  if not attribute in GFF_df.columns]

            # remove missing extra gff attributes
            if any(attributes_absent):

                [ print(f'Warning: Attribute "{attribute}" absent in {GFF_PATH}. Removing.')
                  for attribute in attributes_absent ]

                GFF_info[key]['extras'] = [ attribute 
                                              for attribute in GFF_info[key]['extras'] 
                                              if not attribute in attributes_absent ]

            attributes_duplicated = [ attribute 
                                      for attribute in GFF_info[key]['extras'] 
                                      if attribute in gff_standard_fields]

            # remove extra gff attributes already in standard criteria
            if any(attributes_duplicated):

                [ print(f'Warning: Attribute "{attribute}" already in standard reporting criteria. Removing.')
                  for attribute in attributes_duplicated ]

                GFF_info[key]['extras'] = [ attribute 
                                            for attribute in GFF_info[key]['extras'] 
                                            if not attribute in gff_standard_fields ]


            # Remove URL character encoding from attribute columns (skipping translation if present as this breaks the decoding
            for attribute in GFF_info[key]['extras']:

                if attribute != 'translation':
                
                    GFF_df[attribute] = GFF_df[attribute].fillna('').astype(str).apply(ul.parse.unquote)  # added fix for None datatype



            # extract stanfard gff annotation fields plus any additional annotaiton fields provided
            gff_columns_addback = GFF_df[ [*gff_standard_fields, *GFF_info[key]['extras']] ].copy()

            # fill empty NaN values
            gff_columns_addback.fillna('', inplace=True)
            
            GFF_info[key]['processed'].extend([ gff_columns_addback, GFF_df])

    # extract processed gff info
    gff_basic, gff_pseudo = [ GFF_info[key]['processed'] 
                              for key in sorted(GFF_info.keys()) ]
    
    gff_columns_addback_basic,  attr_to_columns_basic  = gff_basic
    gff_columns_addback_pseudo, attr_to_columns_pseudo = gff_pseudo






    # check seqid consistancy
    
    # read in alignment file object
    pysam_fo = pysam.AlignmentFile( BAM_PATH )

    # extract bam header SQ lines
    bam_header_SQ = pysam_fo.header['SQ']

    # extract SN enteries from each SQ line (conform to string if integer)
    bam_SNs = [ str(name['SN'])
                for name in bam_header_SQ ]

    bam_SNs.sort()

    # extract gff seq ids (conform to string if integer)
    gff_SNs = [ str(name) 
                for name in gff_columns_addback.seq_id.unique() ]

    gff_SNs.sort()

    ##**
    # check if gff and bam sequence IDs are the same or have an overlap
    # if the same run, if an overlap stop and ask if GFF_FORCE should be used of if GFF_FORCE used continue
    # if no overlap stop

    SNs_intersection = set(bam_SNs).intersection(set(gff_SNs))

    SNs_equal = set(bam_SNs) == set(gff_SNs)

    if SNs_equal:

        print('\nGFF and mapped reference sequence IDs match ... proceeding\n')

    else:

        if SNs_intersection:

            print('\nWARNING: GFF & mapping reference sequence IDs are inconsistent')
            print('\nGFF:')
            print(*gff_SNs, sep=", ")
            print('\nBam/Fasta:')
            print(*bam_SNs, sep=", ")
            print('\nOnly the following alignment sequence IDs are present both in the GFF and alignment supplied:')
            print(*sorted(SNs_intersection), sep=", ")
            # [ print(f'-{name}')
            #   for name in sorted(SNs_intersection) ]

            if GFF_FORCE:

                print('\n\t** This sequence ID mismatch has been overridden by --gff_force ... proceeding\n')

            else:
                print('\nERROR: Exiting due to sequence ID mismatch\n')
                print('NOTE: If this sequence ID mismatch is benign e.g. extra/missing minor contigs in the GFF/Fasta, override by using the --gff_force flag\n')
                sys.exit(0)

        else:
            print('\nERROR: Exiting due to no matching sequence IDs')
            print('\nGFF:')
            print(*gff_SNs, sep=", ")
            print('\nBam/Fasta:')
            print(*bam_SNs, sep=", ")
            print('\nNote: Unable to resolve sequence ID problem between the GFF and alignment - check you are using corresponding files\n')
            sys.exit(0)


        # # gff seq ids absent from bam header
        # if not SNs_represented:
        #
        #     # extract missing gff seq ids
        #     missing_SNs = set(gff_SNs).difference(set(bam_SNs))
        #
        #     print('\nWARNING: The following GFF sequence IDs are absent from the alignment file supplied:\n')
        #
        #     [ print(f'-{name}')
        #       for name in sorted(missing_SNs) ]
        #
        #     # gff force flag supplied
        #     if GFF_FORCE:
        #
        #         print('\nSequence ID mismatch will be overridden by --gff_force\n')
        #
        #     else:
        #         print('EXITING: print some message here')
        #         sys.exit(0)

    ##### I EXPECT GFF FORCE WILL RESULT IN PROBLEMS IF GFF SEQ ID MISSING FROM BAM...
    ##** end




    # SAM PROCESSING

    # read in alignment file object
    pysam_fo = pysam.AlignmentFile( BAM_PATH )

    # specify strand ids
    strand = { 0: '+'
              ,1: '-' }
    
    coordinate_delimiter = '\t'

    print('# WRITING COORDINATES TO FILES...')
    with open(bed_file,'w') as f, open(coordinate_file, 'w') as f2:

        # specify coordinate output headers
        coordinate_headers = [ 'ref_name'
                              ,'coord'
                              ,'strand'
                              ,'read_name'
                              ,'read_grp'
                              ,'read_comment' ]
        
        # write headers to read coordinates file
        print( *coordinate_headers,
                sep  = coordinate_delimiter, 
                file = f2 )

        # cycle through alignments
        for read in pysam_fo.fetch():

            # skip unmapped reads
            if read.is_unmapped:
                continue

            # process mapped reads
            else:


                # BED COORDINATES

                # assign appropriate strand id 
                read_strand = strand[int(read.is_reverse)]

                # specify bed coordinates
                read_bed = [ read.reference_name,
                             read.pos,
                             read.reference_end,
                             '.',
                             read.mapping_quality,
                             read_strand,
                             f'# {read.query_name}' ]  # read group added

                # write bed coordinates to output
                print( *read_bed, 
                        sep  = coordinate_delimiter, 
                        file = f )


                # READ COORDINATES

                # extract NM edit distance tag (i.e. number of differences) from alignment fields
                ##### ACCORDING TO MANUAL "If the NM tag is not present, this field will always be 0."
                ##@@ nm_value = read.get_tag('NM')
                ##### HOW DOES THIS AFFECT CALCULATIONS? IS IT FILTERING DATA TOO AGGRESSIVELY?

                ##@@ comment FRACTION MISMATCH check
                # skip if fraction mismatch criteria supplied & reads have excessive mismatches
                # if FRACTION_MISMATCH and (read.query_alignment_length * FRACTION_MISMATCH) > nm_value:
                #     continue
                # ##### originally 0 element of FRACTION_MISMATCH ( i suspect wrong - surely just the entire FRACTION_MISMATCH value - list already unpacked / float non-subscriptable )
                # ##### also surely other way around i.e. nm_value > (read.query_alignment_length * FRACTION_MISMATCH)
                # else:
                ##@@ remove FRACTION MISMATCH check and make else code default
                # specify relevant read coordinate according to strand
                # read_coord = (read.reference_start + 4) if read_strand == '+' else (read.reference_end - 4)
                ##@@ addition of target dup offset parameter default == 4
                read_coord = (read.reference_start + TARGET_DUP_OFFSET) if read_strand == '+' else (read.reference_end - TARGET_DUP_OFFSET)

                # specify read coordinate info
                read_coords = [ read.reference_name,
                                read_coord,
                                read_strand,
                                f'# {read.query_name}',
                                ':'.join(read.query_name.split(':', 4)[:4]),
                                read.get_tag("CO").split(":")[-1] ]  # add fq comment sample id number

                # write read coordinate info to output
                print( *read_coords,
                        sep  = coordinate_delimiter,
                        file = f2 )

    print('# WRITING COORDINATES COMPLETED')
    pysam_fo.close()





    # COORDINATE TO FEATURES

    # read coordinate file in as pandas dataframe
    COORDINATES_df = pd.read_csv( coordinate_file,
                                  sep   = coordinate_delimiter,
                                  dtype = { 'ref_name'     : 'str',
                                            'coord'        : 'int64',
                                            'read_comment' : 'str',
                                            'read_grp'     : 'str' } )

    # bin reads by coordinate & count bin sizes
    COORDINATES_BINNED_df = COORDINATES_df.groupby( by = ['ref_name', 'coord'] ).size().reset_index( name = 'counts' )

    # calculate total insertion sites & total reads
    number_of_insertion_sites = len(COORDINATES_BINNED_df)
    number_of_reads_mapped    = COORDINATES_BINNED_df['counts'].sum()

    # calculate minimum / maximum / median / mean read support at insertion sites
    min_reads_at_site         = COORDINATES_BINNED_df['counts'].min()
    max_reads_at_site         = COORDINATES_BINNED_df['counts'].max()
    median_reads_at_site      = round(COORDINATES_BINNED_df['counts'].median(), 2)
    mean_insertion_site_depth = round(number_of_reads_mapped / number_of_insertion_sites, 2)

    # remove coordinates below specified depth
    COORDINATES_BINNED_df.query( f' counts >= {MIN_DEPTH} ', inplace=True )

    FILTERED_df = COORDINATES_BINNED_df.copy()

    # insert additional data columns
    FILTERED_df.eval(
        f''' 
        source       = "pimms2"
        feature_type = "misc_feature"
        strand       = "."
        phase        = "."
        score        =  counts
        start        =  coord
        stop         =  coord
        info         = "note=insertion;"
        ''', inplace=True)

    # specify gff coordinates path
    gff_coordinates = os.path.join( db_rdir, 
                                    f'{RESULTS_LABEL}_pimms_insert_coordinates.gff')

    # specify relevant columns
    relevant_columns = [ 'ref_name',
                         'source',
                         'feature_type',
                         'start',
                         'stop',
                         'score',
                         'strand',
                         'phase',
                         'info' ]

    # write gff coordinates to output
    FILTERED_df[ [ *relevant_columns ] ].to_csv( gff_coordinates, 
                                                 index  = False,
                                                 sep    = coordinate_delimiter, 
                                                 header = False )

    # calculate gap size between insert sites within same ref_name
    COORDINATES_BINNED_df['between_insertion_gap'] = COORDINATES_BINNED_df.groupby( by='ref_name' )['coord'].diff()

    # calculate minimum / maximum / median / mean gap size between insertion sites
    min_between_insertion_gap    = COORDINATES_BINNED_df['between_insertion_gap'].min()
    max_between_insertion_gap    = COORDINATES_BINNED_df['between_insertion_gap'].max()
    median_between_insertion_gap = COORDINATES_BINNED_df['between_insertion_gap'].median()
    mean_between_insertion_gap   = round(COORDINATES_BINNED_df['between_insertion_gap'].mean(), 2)
    
    # specify sql to extract standard gff attributes
    sqlcode = '''
        select COORDINATES_BINNED_df.*
        ,attr_to_columns_basic.*
        from attr_to_columns_basic
        left join COORDINATES_BINNED_df
        on COORDINATES_BINNED_df.coord between attr_to_columns_basic.start and attr_to_columns_basic.end
        where COORDINATES_BINNED_df.ref_name like '%' || attr_to_columns_basic.seq_id || '%'
        '''
    # N.B. wierd sqlite concatenation + >> || ,  '%' == wildcard double check effect of this
    # this line should allow multi contig files

    # extract relevant gff features near insertion sites
    COORDINATES_GFF_df = ps.sqldf(sqlcode, locals())
    
    # calculate forward strand position percentile
    posn_as_percentile_neg = ( ( COORDINATES_GFF_df.coord - COORDINATES_GFF_df.start ) + 1 ) / ( COORDINATES_GFF_df.feat_length / 100 )

    # calculate reverse strand position percentile
    posn_as_percentile_pos = ( ( COORDINATES_GFF_df.end   - COORDINATES_GFF_df.coord ) + 1 ) / ( COORDINATES_GFF_df.feat_length / 100 )

    # select relevant positions percentile according to strand orientation
    COORDINATES_GFF_df['posn_as_percentile'] = posn_as_percentile_neg.where( cond  = COORDINATES_GFF_df.strand == '+',
                                                                             other = posn_as_percentile_pos ).round( decimals = 1 )

    # Important fix group by doesn't work -- any rows with nan values get dropped *yikes very bad!!!!*
    COORDINATES_GFF_df.fillna( value   = '', 
                            inplace = True )

    # groupby gff fields (filtered to remove absent & duplicate extra fields) & calculate summary metrics
    pimms_result_table = COORDINATES_GFF_df.groupby( [ *gff_columns_addback_basic.columns ] ).agg( num_insertions_mapped_per_feat  = ( 'counts',             'sum'   ),
                                                                                                   num_insert_sites_per_feat       = ( 'counts',             'count' ),
                                                                                                   first_insert_posn_as_percentile = ( 'posn_as_percentile', 'min'   ),
                                                                                                   last_insert_posn_as_percentile  = ( 'posn_as_percentile', 'max'   ) ).reset_index()

    # calculate additional summary metrics
    pimms_result_table.eval(
        f'''
         num_insert_sites_per_feat_per_kb = ( num_insert_sites_per_feat      /   feat_length ) * 1000
         NRM_score                        = ( num_insertions_mapped_per_feat / ( feat_length   / 1000 ) ) / ( {number_of_reads_mapped} / 1e6 )
         NIM_score                        = ( num_insert_sites_per_feat      / ( feat_length   / 1000 ) ) / ( {number_of_reads_mapped} / 1e6 )
         ''', inplace=True)

    pimms_result_table = pimms_result_table.round({ 'num_insert_sites_per_feat_per_kb' : 2,
                                                    'NRM_score'                        : 2,
                                                    'NIM_score'                        : 2 })

    # extract metric columns
    metric_columns = set(pimms_result_table.columns).difference(set(gff_columns_addback_basic.columns))

    # get insertion metric column names
    if RESULTS_LABEL:
        
        # add prefix to metric columns
        renamed_columns = { column: f'{RESULTS_LABEL}_{column}' 
                            for column in metric_columns }

        # rename metric columns
        pimms_result_table.rename( columns = renamed_columns,
                                   inplace = True)
    
        # update metric columns
        metric_columns = [ column 
                           for column in renamed_columns.values() ]
    
    # merge insertion metrics into previous dataframe
    pimms_result_table_full = pd.merge( left  = gff_columns_addback_basic, 
                                        right = pimms_result_table, 
                                        how   = 'left'  )

    # fill empty values where insertion metrics absent
    pimms_result_table_full.fillna( value   = { metric: int(0)
                                                for metric in metric_columns }, 
                                    inplace = True )



    # COORDINATES -> FEATURES (REPS)

    if SKIP_REPS:
        
        print("not processing read groups as replicates\n")

    else:
        
        print("processing read groups as replicates for illumina data\n")

        # extract unique read groups from coordinates file
        read_grps = sorted(COORDINATES_df.read_grp.unique())
        print(f'\n{len(read_grps)} readgroups found:')
        [ print(f'-{group}') for group in read_grps ]

        # extract unique read comments from coordinates file
        read_comments = sorted(COORDINATES_df.read_comment.unique())

        ##**
        # fix to override reps if many comments found due to poor parsing of old solexa style fastq header lines
        if (len(read_comments) > 12):
            print("not processing read groups as replicates\n")
            print("Warning: Unable to resolve different samples in fastq/sam/bam data (apparently too many?), continuing without replicates")
            print("Note: This may be due to old style Illumina header lines")
            pool_comment_criteria = False
            pool_group_criteria = False

        else:
            print(f'\n{len(read_comments)} sample comments found:')
            [ print(f'-{comment}') for comment in read_comments ]

            # specify respective pooling criteria
            pool_comment_criteria = (len(read_grps) < len(read_comments)) & (len(read_comments) >= 3)

            pool_group_criteria   = (len(read_grps) >= len(read_comments)) & (len(read_grps)     >= 3) ##**  # change to capture if number of read_grps amd read_comments is same
        ##** end

        # skip replicates        
        if not pool_comment_criteria and not pool_group_criteria:

            print('\nWarning: Unable to resolve >= 3 samples in fastq/sam/bam data, continuing without replicate insertion counts')
            print('N.B: If this is an error the software may need updating to recognise novel fastq naming conventions')
        
        # pool replicates
        else:

            # set read comments for pooling
            if pool_comment_criteria:
                
                field_drop   = 'read_grp'
                field_pool   = 'read_comment'
                field_unique = read_comments


            # set read groups for pooling
            elif pool_group_criteria:

                field_drop   = 'read_comment'
                field_pool   = 'read_grp'
                field_unique = read_grps

            # drop irrelevant column
            COORDINATES_df.drop( labels  = field_drop, 
                                 axis    = 1,
                                 inplace = True)

            # rename column for pooling
            COORDINATES_df.rename( columns = { field_pool: 'sample_info' },
                                   inplace = True )

            # calculate pool sizes
            COORDINATES_df = COORDINATES_df.groupby( by = ['ref_name', 'coord', 'sample_info'] ).size().reset_index( name = RESULTS_LABEL )
            
            # report pool counts
            ##### REALLY!?!?! SURELY MORE BECAUSE INLUCDING REF_NAME & COORD IN GROUPBY (I.E. MORE UNIQUE COMBINATIONS OF REF_NAME + COORD + SAMPLE_INFO)
            print(f'{len(field_unique)} sample replicates/mutant pools established')
            ##### SUBSTITUTE WITH DATAFRAME LENGTH?
            
            # reshape dataframe with pool size columns for each unique sample_info
            COORDINATES_df = COORDINATES_df.copy( deep=False ).pivot_table( index      = [ 'ref_name', 'coord' ],
                                                                            columns    = [ 'sample_info' ],
                                                                            values     = [ RESULTS_LABEL ],
                                                                            fill_value = 0 ).reset_index()
            
            # replace sample_info names with numerical labels
            COORDINATES_df.rename( columns = { name : f'MP{i}' 
                                               for i, name in enumerate(field_unique,1) }, 
                                   inplace = True )
            
            # flatten column indexes with hiererchial label included as prefix
            COORDINATES_df.columns = [ '_'.join(col).rstrip('_ ')
                                             for col in COORDINATES_df.columns.to_flat_index() ]

            sqlcode = '''
                select COORDINATES_df.*
                ,attr_to_columns_basic.seq_id
                ,attr_to_columns_basic.start
                ,attr_to_columns_basic.end
                from attr_to_columns_basic
                left join COORDINATES_df
                on COORDINATES_df.coord between attr_to_columns_basic.start and attr_to_columns_basic.end
                where COORDINATES_df.ref_name like '%' || attr_to_columns_basic.seq_id || '%'
                '''

            # wierd sqlite concatenation + >> || ,  '%' == wildcard double check effect of this
            # this line should allow multi contig files

            # extract relevant gff features near insertion sites
            COORDINATES_GFF_df = ps.sqldf(sqlcode, locals())

            # remove ref_name & coord columns
            COORDINATES_GFF_df.drop(labels  = ['ref_name','coord'],
                                    axis    = 1, 
                                    inplace = True)

            # groupby by gff feature coordinates & sum
            COORDINATES_GFF_df = COORDINATES_GFF_df.groupby(['seq_id', 'start', 'end']).sum().reset_index()

            # merge feature metrics into previous dataframe
            pimms_result_table_full = pd.merge( left  = pimms_result_table_full, 
                                                right = COORDINATES_GFF_df, 
                                                on  = ["seq_id", "start", "end"],
                                                how   = 'left' )

            # fill empty values where feature metrics absent
            pimms_result_table_full.fillna( value   = { metric: int(0)
                                                        for metric in COORDINATES_GFF_df.columns
                                                        if metric.startswith(RESULTS_LABEL) }, 
                                            inplace = True )


    # pseudogene features present
    if not gff_columns_addback_pseudo.empty:
    

        # extract pseudogene features
        pseudo_entries = pimms_result_table_full.locus_tag.isin(gff_columns_addback_pseudo['locus_tag'])

        # tag feature type to indicuate pseudogene
        pimms_result_table_full.loc[pseudo_entries, "type"] = pimms_result_table_full['type'] + '_pseudo'


    # detect output mode from provided extension
    output_modes = [ OUT_FORMAT == extension 
                    for extension in ('xlsx', 'csv', 'tsv') ]

    xlsx_mode, csv_mode, tsv_mode = output_modes

    # write countinfo outputs

    if xlsx_mode:

        writer = pd.ExcelWriter( countinfo_file, 
                                 engine='xlsxwriter' )
        
        # Convert the dataframe to an XlsxWriter Excel object.
        pimms_result_table_full.to_excel( writer, 
                                          sheet_name='PIMMS2_result', 
                                          index=False )

        # Close the Pandas Excel writer and output the Excel file.
        writer.save()

    if csv_mode:

        pimms_result_table_full.to_csv( countinfo_file, 
                                        index=False, 
                                        sep=',' )

    if tsv_mode:

        pimms_result_table_full.to_csv( countinfo_file, 
                                        index=False, 
                                        sep='\t' )










#######################
# TABLE_MERGE COMMAND #
#######################

if MODE_SELECTED == TABLE_MERGE:
    
    # filter absent input tables
    input_tables = [ inputs
                     for inputs in [ parsed_args[0].xlsx
                                    ,parsed_args[0].csv
                                    ,parsed_args[0].tsv ]
                     if inputs ]

    # decide whether to proceed with merging tables
    if not input_tables:
        print("\nERROR: Result tables absent for merging\n")
    else:

        # unpack present input tables
        (table1, table2), = input_tables, = input_tables

        print(f'Merging:')
        [ print(f'\t-{table}') for table in input_tables ]

        # detect merge mode from table extension
        merge_modes = [ all( table.endswith(extension) 
                             for table in input_tables ) 
                        for extension in ('xlsx', 'csv', 'tsv') ]

        xlsx_mode, csv_mode, tsv_mode = merge_modes


    # READ INPUTS

        # convert to pandas dataframe
        if xlsx_mode:

            result_dfs = [ pd.read_excel(table, engine="openpyxl")
                           for table in input_tables ]
            
        if csv_mode:

            result_dfs = [ pd.read_csv(table).replace('"', '', regex=True)
                           for table in input_tables ]

        if tsv_mode:

            result_dfs = [ pd.read_csv(table, sep="\t")
                           for table in input_tables ]


    # MERGE TABLES

        # extract individual dataframes
        df1, df2 = result_dfs

        # merge pandas dataframes 
        merged_dfs = pd.merge( left=df1
                                       , right=df2 )
            

    # WRITE OUTPUTS

        if xlsx_mode:

            writer = pd.ExcelWriter('merged_result.xlsx', engine='xlsxwriter')
            # Convert the dataframe to an XlsxWriter Excel object.
            merged_dfs.to_excel(writer, sheet_name='PIMMS2_merged_result', index=False)
            # Close the Pandas Excel writer and output the Excel file.
            writer.save()
            
        if csv_mode:

            merged_dfs.to_csv('merged_result.csv', index=False, sep=',')
        
        if tsv_mode:
            
            merged_dfs.to_csv('merged_result.txt', index=False, sep="\t")
