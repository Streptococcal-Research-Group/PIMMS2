import os
import sys
import datetime
import multiprocessing
import pysam
import pandas as pd
import gffpandas.gffpandas as gffpd
import pandasql as ps
import pysam
import urllib as ul
import configargparse


# from configargparse import _SubParsersAction
# import configparser
## pandas + 'xlsxwriter'

# import configparser
#
# sub get_position_as_percentile_pos{
# my $in = $_;
# my $in_locus  = $locus_lookup{$_};
# my $in_length = $length_lookup{$in_locus};
# my $in_start = $start_lookup{$in_locus};
# my $percentile =  sprintf("%.1f", ((($in-$in_start)+1) / ($in_length/100)));
# return $percentile;
# }
#
# sub get_position_as_percentile_neg{
# my $in = $_;
# my $in_locus  = $locus_lookup{$_};
# my $in_length = $length_lookup{$in_locus};
# my $in_stop = $stop_lookup{$in_locus};
# my $percentile =  sprintf("%.1f", ((($in_stop-$in)+1) / ($in_length/100)));
# return $percentile;
# }
class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end


# main function
# def main():
from gffpandas.gffpandas import Gff3DataFrame

ap = configargparse.ArgumentParser(  # description='PIMMS2 sam/bam processing',
    prog="demo_pimms2_sam_extract",
    add_config_file_help=False,
    config_file_parser_class=configargparse.DefaultConfigFileParser,
    epilog="\n\n*** N.B. This is a development version ***\n \n ",
    description='''this description
                                   was indented weird
                                   but that is okay'''
)
ap.add_argument('--version', action='version', version='%(prog)s 2.0.1 demo')
modes = ap.add_subparsers(parser_class=configargparse.ArgParser)

# modes.required = False
findflank = modes.add_parser("find_flank", help="Mode: filter fastq files to find insertion site flanking sequence")
samcoords = modes.add_parser("sam_extract", add_config_file_help=False,
                             help="Mode: extract insertion site coordinates from sam file",
                             description="Args that start with '--' (eg. --sam) can also be set in a config file (specified via -c)")
otherstuff = modes.add_parser("other_stuff", help='Mode: do other good PIMMS2 related stuff')

samcoords.add_argument("-c", "--config", required=False, is_config_file=True,  # dest='config_file',
                       # type=str, default='',
                       metavar='pimms2.config',
                       help="use parameters from config file")
samcoords.add_argument("--sam", required=True, nargs=1, metavar='pimms.sam/bam',
                       help="sam/bam file of mapped IS flanking sequences ")
samcoords.add_argument("--label", required=False, nargs=1, metavar='condition_name', default=[''],
                       help="text tag to add to results file")
samcoords.add_argument("--mismatch", required=False, nargs=1, type=float, metavar='float', default=[None],
                       choices=[Range(0.0, 0.2)],
                       help="fraction of permitted mismatches in mapped read ( 0 <= float < 0.2 [no filter]")
samcoords.add_argument("--min_depth", required=False, nargs=1, type=int, default=1, metavar='int',
                       help="minimum read depth at insertion site")
samcoords.add_argument("--gff", required=True, nargs=1, type=str, default='', metavar='genome.gff',
                       help="GFF3 formatted file to use\n(note fasta sequence present in the file must be deleted before use)")
samcoords.add_argument("--gff_extra", required=False, nargs=1, type=str, default='', metavar="'x,y,z'",
                       help="comma separated list of extra fields to include from the GFF3 annotation\ne.g. 'ID,translation,note' ")

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

# print((vars(parsed_args)))
# exit()
# do config parsing
# config_file = parsed_args.config_file[0]

# construct config parser
# p2config = configparser.ConfigParser()

# p2config.read(config_file)
# print("extra gff fields: " + (parsed_args.gff_extra))
if parsed_args[0].gff_extra:
    # strip any formatting quotes and turn comma separated string into a list of fields
    gff_extra = parsed_args[0].gff_extra[0].strip("'\"").split(',')
else:
    gff_extra = []

print("extra gff fields: " + str(gff_extra))

# exit()
use_fraction_mismatch = True
min_depth_cutoff = parsed_args[0].min_depth[0]
# permitted_mismatch = 6
fraction_mismatch = parsed_args[0].mismatch[0]
# sam_file = "0140J_test.IN.PIMMS.RX.rmIS.sub0_min25_max50.ngm_sensitive.90pc.bam"

# sam_file = "UK15_Media_RX_pimms2out_trim60_v_pGh9_UK15.bam"
# sam_file = "UK15_Blood_RX_pimms2out_trim60_v_pGh9_UK15.bam"
sam_file = parsed_args[0].sam[0]
sam_stem, sam_ext = os.path.splitext(os.path.basename(sam_file))
sam_stem = sam_stem + '_md' + str(min_depth_cutoff) + '_mm' + str(fraction_mismatch or '')
# gff_file = '/Users/svzaw/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/S.uberis_0140J.gff3'
# gff_file = 'UK15_assembly_fix03.gff'
gff_file = parsed_args[0].gff[0]
gff_feat_type = ['CDS', 'tRNA', 'rRNA']

# exit()

annotation = gffpd.read_gff3(gff_file)
annotation = annotation.filter_feature_of_type(gff_feat_type)
# break 9th gff column key=value pairs down to make additional columns
attr_to_columns = annotation.attributes_to_columns()
# attr_to_columns = attr_to_columns[attr_to_columns['type'] == gff_feat_type]  # filter to CDS
# attr_to_columns = attr_to_columns.filter_feature_of_type(gff_feat_type)  # filter to gff_feat_type = ['CDS', 'tRNA', 'rRNA']
# add feature length column and do some cleanup
attr_to_columns = attr_to_columns.assign(
    feat_length=(attr_to_columns.end - attr_to_columns.start + 1)).dropna(axis=1,
                                                                          how='all').drop(columns=['attributes'])

# remove RFC 3986 % encoding from product (gff3 attribute)
# attr_to_columns = attr_to_columns.assign(product_nopc=attr_to_columns['product'].apply(ul.parse.unquote)).drop(
#    columns=['product']).rename(columns={'product_nopc': 'product'})

attr_to_columns['product'] = attr_to_columns['product'].apply(ul.parse.unquote)

for field in gff_extra:
    attr_to_columns[field] = attr_to_columns[field].apply(ul.parse.unquote)

gff_columns_addback = attr_to_columns[['seq_id',
                                       'locus_tag',
                                       'gene',
                                       'start',
                                       'end',
                                       'feat_length',
                                       'product'] + gff_extra]  # add extra fields from gff
# attr_to_columns = attr_to_columns.dropna(axis=1, how='all')
samfile = pysam.AlignmentFile(sam_file)  # without , "rb" should auto detect sam or bams
open(sam_stem + ".bed", 'w').close()
f = open(sam_stem + ".bed", "a")

STRAND = ["+", "-"]
for read in samfile.fetch():
    STR = STRAND[int(read.is_reverse)]
    BED = [read.reference_name, read.pos, read.reference_end, ".", read.mapping_quality, STR, '# ' + read.query_name]
    f.write('\t'.join([str(i) for i in BED]))
    f.write('\n')
f.close()

open(sam_stem + ".insert_coords.txt", 'w').close()
f2 = open(sam_stem + ".insert_coords.txt", "a")
f2.write('\t'.join([str(i) for i in ['ref_name', 'coord', 'strand', 'read_name']]))
f2.write('\n')
STRAND = ["+", "-"]
for read in samfile.fetch():
    if read.is_unmapped:
        continue
    #   if read.infer_query_length == 25:
    #   print(str(read.get_tag('NM')))
    NM_value = read.get_tag('NM')
    # if not use_fraction_mismatch and NM_value > permitted_mismatch:
    #    continue
    if fraction_mismatch:  # and NM_value > 0:
        if (read.query_alignment_length * fraction_mismatch[0]) > NM_value:
            continue
    STR = STRAND[int(read.is_reverse)]  # coverts is_reverse boolean into + or - strings
    if STR == '+':
        COORDS = [read.reference_name, (read.reference_start + 4), STR, '# ' + read.query_name]
        f2.write('\t'.join([str(i) for i in COORDS]))
        f2.write('\n')
    if STR == '-':
        COORDS = [read.reference_name, (read.reference_end - 4), STR, '# ' + read.query_name]
        f2.write('\t'.join([str(i) for i in COORDS]))
        f2.write('\n')
f2.close()

coord_df = pd.read_csv(sam_stem + ".insert_coords.txt", sep='\t', dtype={'ref_name': "str", 'coord': "int64"})
coord_counts_df = coord_df.groupby(['ref_name', 'coord']).size().reset_index(name='counts')

number_of_insertion_sites = len(coord_counts_df)
number_of_reads_mapped = coord_counts_df['counts'].sum()
min_reads_at_site = coord_counts_df['counts'].min()
max_reads_at_site = coord_counts_df['counts'].max()
median_reads_at_site = round(coord_counts_df['counts'].median(), 2)

mean_insertion_site_depth = round(number_of_reads_mapped / number_of_insertion_sites, 2)

coord_counts_df = coord_counts_df[coord_counts_df['counts'] > min_depth_cutoff]

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
## wierd sqlite concatenation + >> || ,  '%' == wildcard double check effect of this
# this line should allow multi contig files

coords_join_gff = ps.sqldf(sqlcode, locals())
coords_join_gff.count
## add position as percentile (needs manual confirmation)
coords_join_gff = coords_join_gff.assign(
    # python/pandas implementation of PIMMS.pl code to derive insert position as percentile of gene length
    # sprintf("%.1f", ((($in-$in_start)+1) / ($in_length/100)));
    # sprintf("%.1f", ((($in_stop-$ in) + 1) / ($in_length / 100)));
    posn_as_percentile=(((coords_join_gff.coord - coords_join_gff.start) + 1) / (
            coords_join_gff.feat_length / 100)).where(
        coords_join_gff.strand == '+', ((coords_join_gff.end - coords_join_gff.coord) + 1) / (
                coords_join_gff.feat_length / 100))).round({"posn_as_percentile": 1})
print(list(attr_to_columns.columns.values))

pimms_result_table = coords_join_gff.groupby(
    ['seq_id', 'locus_tag', 'gene', 'start', 'end', 'feat_length', 'product'] + gff_extra).agg(
    num_reads_mapped_per_feat=('counts', 'sum'),
    num_insert_sites_per_feat=('counts', 'count'),
    first_insert_posn_as_percentile=('posn_as_percentile', 'min'),
    last_insert_posn_as_percentile=('posn_as_percentile', 'max')
).reset_index()

pimms_result_table = pimms_result_table.assign(num_insert_sites_per_feat_per_kb=(
        (pimms_result_table.num_insert_sites_per_feat / pimms_result_table.feat_length) * 1000),
    NRM_score=((pimms_result_table.num_reads_mapped_per_feat / (
            pimms_result_table.feat_length / 1000)) / (
                       number_of_reads_mapped / 1e6)),
    NIM_score=((pimms_result_table.num_insert_sites_per_feat / (
            pimms_result_table.feat_length / 1000)) / (
                       number_of_reads_mapped / 1e6))
).round({'num_insert_sites_per_feat_per_kb': 2, 'NRM_score': 2, 'NIM_score': 2})
print(list(pimms_result_table.columns.values))

pimms_result_table = pimms_result_table[['seq_id',
                                         'locus_tag',
                                         'gene',
                                         'start',
                                         'end',
                                         'feat_length',
                                         'product'] + gff_extra +
                                        ['num_reads_mapped_per_feat',
                                         'num_insert_sites_per_feat',
                                         'num_insert_sites_per_feat_per_kb',
                                         'first_insert_posn_as_percentile',
                                         'last_insert_posn_as_percentile',
                                         'NRM_score',
                                         'NIM_score']]

print(list(pimms_result_table.columns.values))

# pimms_result_table_full gff_columns_addback

NAvalues = {'num_reads_mapped_per_feat': int(0),
            'num_insert_sites_per_feat': int(0),
            'num_insert_sites_per_feat_per_kb': int(0),
            'first_insert_posn_as_percentile': int(0),
            'last_insert_posn_as_percentile': int(0),
            'NRM_score': int(0),
            'NIM_score': int(0)}
pimms_result_table_full = pd.merge(gff_columns_addback, pimms_result_table, how='left').fillna(value=NAvalues)

# if set add prefix to columns
if parsed_args[0].label[0]:
    label_cols = pimms_result_table_full.columns[pimms_result_table_full.columns.isin(['num_reads_mapped_per_feat',
                                                                                       'num_insert_sites_per_feat',
                                                                                       'num_insert_sites_per_feat_per_kb',
                                                                                       'first_insert_posn_as_percentile',
                                                                                       'last_insert_posn_as_percentile',
                                                                                       'NRM_score',
                                                                                       'NIM_score'])]
    pimms_result_table_full.rename(columns=dict(zip(label_cols, parsed_args[0].label[0] + '_' + label_cols)),
                                   inplace=True)

pimms_result_table_full.to_csv(sam_stem + "_countinfo_tab.txt", index=False, sep='\t')
pimms_result_table_full.to_csv(sam_stem + "_countinfo_tab.csv", index=False, sep=',')

writer = pd.ExcelWriter(sam_stem + '_countinfo.xlsx', engine='xlsxwriter')

# Convert the dataframe to an XlsxWriter Excel object.
pimms_result_table_full.to_excel(writer, sheet_name='PIMMS2_result', index=False)

# Close the Pandas Excel writer and output the Excel file.
writer.save()

# num_insert_sites_per_feat_per_kb=('counts', '(count / feat_length) *100')

# ['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'Dbxref', 'ID', 'Name', 'Note', 'Parent', 'gbkey', 'gene', 'locus_tag', 'product', 'protein_id', 'pseudo', 'transl_table', 'feat_length']
# 'seq_id', 'locus_tag', 'gene', 'start', 'end', 'feat_length', 'product',
# Locus	Locus information from GTF file as determined by -g
# Gene	Gene information from GTF file as determined by -g
# Start	Gene start position from GTF file as determined by -g
# Stop	Gene stop position from GTF file as determined by -g
# CDS length	Locus length
# Product	Product information from GTF file as determined by -g
# Number of mutations	Total number of insertions mapped within this locus
# Number of unique mutations (>= [-m] x coverage)	Total number of unique insertion positions mapped within this locus
# Unique mutations per1KbCDS (>= [-m] x coverage)	Total number of unique insertion positions mapped within this locus per kb of each locus
# First unique insertion centile position (>= [-m] x coverage)	For each unique insertion position, the centile position of that locus is determined. The first (lowest centile) position is returned.
# Last unique insertion centile position(>= [-m] x coverage)	For each unique insertion position, the centile position of that locus is determined. The last (greatest centile) position is returned.
# Normalised Reads Mapped (NRM score)	(total number of reads mapped per gene / length of gene in KB) / (total mapped read count / 10 ^ 6)
# Normalised Insertions Mapped (NIM score)	(total number of unique insertions mapped per gene/length of gene in KB)/(total mapped read count/10^6)

#    if __name__ == '__main__':
#       main()
