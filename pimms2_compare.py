import os
import sys
import datetime
# import multiprocessing
# import pysam
import pandas as pd
# import gffpandas.gffpandas as gffpd
import pandasql as ps

import configargparse

ap = configargparse.ArgumentParser(  # description='PIMMS2 sam/bam processing',
    prog="demo_pimms2_compare",
    add_config_file_help=False,
    config_file_parser_class=configargparse.DefaultConfigFileParser,
    epilog="\n\n*** N.B. This is a development version ***\n \n ",
    description='''description here'''
)
ap.add_argument('--version', action='version', version='%(prog)s 2.0.1 demo')
modes = ap.add_subparsers(parser_class=configargparse.ArgParser)

# modes.required = False
findflank = modes.add_parser("find_flank", add_config_file_help=False,
                             help="Mode: filter fastq files to find insertion site flanking sequence")
samcoords = modes.add_parser("sam_extract", add_config_file_help=False,
                             help="Mode: extract insertion site coordinates from sam file",
                             description="Args that start with '--' (eg. --sam) can also be set in a config file (specified via -c)")
compare = modes.add_parser("compare", add_config_file_help=False,
                           help='Mode: do other good PIMMS2 related stuff')

compare.add_argument("-c", "--config", required=False, is_config_file=True,  # dest='config_file',
                     # type=str, default='',
                     metavar='pimms2.config',
                     help="use parameters from config file")
compare.add_argument("--result1", required=True, nargs=1, metavar='pimms.sam/bam',
                     help="pimms result files condition 1 ")
compare.add_argument("--result2", required=True, nargs=1, metavar='pimms.sam/bam',
                     help="pimms result files condition 2 ")
compare.add_argument("--label1", required=True, nargs=1, metavar='condition_one',
                     help="text tag for condition one")
compare.add_argument("--label2", required=True, nargs=1, metavar='condition_two',
                     help="text tag for condition two")

# parsed_args = ap.parse_args()
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

= "UK15_Media_RX_pimms2out_trim60_v_pGh9_UK15_md1_mm_countinfo_tab.txt"
result2 = "UK15_Blood_RX_pimms2out_trim60_v_pGh9_UK15_md1_mm_countinfo_tab.txt"
label1 = "media"
label2 = "blood"

result1_df = pd.read_csv(result1, sep='\t')
result2_df = pd.read_csv(result2, sep='\t')
result1_df = result1_df.filter(regex='^media|seq_id|locus_tag', axis=1)
# result1_df = result1_df[result1_df[label1 + '_num_reads_mapped_per_feat'] > 0]
# result2_df = result2_df[result2_df[label2 + '_num_reads_mapped_per_feat'] > 0]

result1_2_df = pd.DataFrame.merge(result1_df, result2_df, how='outer', on=["seq_id", "locus_tag"])
result1_2_filt_df = result1_2_df[(result1_2_df[label1 + '_' + 'num_reads_mapped_per_feat'] > 0) |
                                 (result1_2_df[label2 + '_' + 'num_reads_mapped_per_feat'] > 0)]

result1_2_filt_df = result1_2_filt_df.insert(2, NIM_diff=(
        (pimms_result_table.num_insert_sites_per_feat / pimms_result_table.feat_length) * 1000
