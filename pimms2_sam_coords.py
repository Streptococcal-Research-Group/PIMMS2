import os
import sys
import re
# import datetime
# import multiprocessing
import pandas as pd
# pd.set_option('mode.chained_assignment', 'raise')

import gffpandas.gffpandas as gffpd
import pandasql as ps
import pysam
import urllib as ul
import configargparse
from pathlib import Path
import shutil


# test commit


# import warnings

# from configargparse import _SubParsersAction
# import configparser
## pandas + 'xlsxwriter'

# import configparser
#

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise configargparse.ArgumentTypeError("{0} does not exist".format(x))
    return x






def process_gff(gff_file, gff_feat_type, gff_extra):
    annotation = gffpd.read_gff3(gff_file)
    annotation = annotation.filter_feature_of_type(gff_feat_type)

    gff_stem, gff_ext = os.path.splitext(os.path.basename(gff_file))

    annotation.to_gff3(gff_stem + '_pimms_features.gff')
    # break 9th gff column key=value pairs down to make additional columns
    attr_to_columns = annotation.attributes_to_columns()

    if attr_to_columns.empty:
        # return empty dataframes if no rows of required type are found print('attr_to_columns is empty!')
        return attr_to_columns, attr_to_columns
    # attr_to_columns = attr_to_columns[attr_to_columns['type'] == gff_feat_type]  # filter to CDS
    # attr_to_columns = attr_to_columns.filter_feature_of_type(gff_feat_type)  # filter to gff_feat_type = ['CDS', 'tRNA', 'rRNA']
    # add feature length column and do some cleanup
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
        attr_to_columns['product'] = attr_to_columns['product'].fillna('').astype(str).apply(
            ul.parse.unquote)  # added fix for None datatype
    ## fix to skip requested extra gff annotation field if not present in GFF
    drop_gff_extra = []
    for field in gff_extra:
        if field not in attr_to_columns:
            print("Warning: Unable to find '" + field + "' in " + str(gff_file) + ' file, continuing...')
            drop_gff_extra.append(field)

    gff_extra = [item for item in gff_extra if item not in drop_gff_extra]

    ## Remove URL character encoding from columns  (skipping translation if present as this breaks the decoding
    for field in gff_extra:
        if field == 'translation':
            continue
        else:
            attr_to_columns[field] = attr_to_columns[field].fillna('').astype(str).apply(
                ul.parse.unquote)  # added fix for None datatype

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
    # gff_columns_addback.fillna(dict.fromkeys(['seq_id',
    #                                           'locus_tag',
    #                                           'type',
    #                                           'gene'], 0), inplace=True)

    data_top = gff_columns_addback.head()
    print(data_top)
    data_top2 = attr_to_columns.head()
    print(data_top2)

    return gff_columns_addback, attr_to_columns
    # end process_gff()


def modify_sam_stem(sam_file):
    sam_stem, sam_ext = os.path.splitext(os.path.basename(sam_file))
    sam_stem: str = sam_stem + '_md' + str(min_depth_cutoff) + '_mm' + str(fraction_mismatch or '')
    return sam_stem


def process_sam(sam_file):
    sam_stem = modify_sam_stem(sam_file)
    samfile = pysam.AlignmentFile(sam_file)  # without , "rb" should auto detect sam or bams

    open(sam_stem + ".bed", 'w').close()
    f = open(sam_stem + ".bed", "a")

    STRAND = ["+", "-"]
    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        STR = STRAND[int(read.is_reverse)]
        BED = [read.reference_name, read.pos, read.reference_end, ".", read.mapping_quality, STR,
               '# ' + read.query_name]  # read group added
        f.write('\t'.join([str(i) for i in BED]))
        f.write('\n')
    f.close()
    print('#BED')
    samfile.close()

    samfile = pysam.AlignmentFile(sam_file)
    open(sam_stem + "_insert_coords.txt", 'w').close()
    f2 = open(sam_stem + "_insert_coords.txt", "a")
    f2.write('\t'.join([str(i) for i in ['ref_name', 'coord', 'strand', 'read_name', 'read_grp']]))
    f2.write('\n')
    STRAND = ["+", "-"]
    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        # print(read.query_name + '.')
        # continue
        NM_value = read.get_tag('NM')
        if fraction_mismatch:  # and NM_value > 0:
            if (read.query_alignment_length * fraction_mismatch[0]) > NM_value:
                continue
        STR = STRAND[int(read.is_reverse)]  # coverts is_reverse boolean into + or - strings
        # print(STR)
        #continue
        if STR == '+':
            COORDS = [read.reference_name, (read.reference_start + 4), STR, '# ' + read.query_name,
                      ':'.join(read.query_name.split(':', 4)[:4])]
            f2.write('\t'.join([str(i) for i in COORDS]))
            f2.write('\n')
        if STR == '-':
            COORDS = [read.reference_name, (read.reference_end - 4), STR, '# ' + read.query_name,
                      ':'.join(read.query_name.split(':', 4)[:4])]
            f2.write('\t'.join([str(i) for i in COORDS]))
            f2.write('\n')
    f2.close()
    samfile.close()
    print('#COORDS')
    # end process_sam()


def seqID_consistancy_check(mygffcolumns, my_sam):
    af = pysam.AlignmentFile(my_sam)
    sam_seq_ID_list = [name['SN'] for name in af.header['SQ']]
    gff_seq_ID_list = mygffcolumns.seq_id.unique().tolist()
    sam_seq_ID_list.sort()
    gff_seq_ID_list.sort()
    if sam_seq_ID_list == gff_seq_ID_list:
        print('GFF & mapping reference sequence IDs match')
    elif parsed_args[0].gff_force:
        print('\nWARNING: GFF & mapping reference sequence IDs are inconsistent. \n' +
              'sequence ID mismatch overridden by --gff_force\ngff:\n' +
              str(gff_seq_ID_list) + '\nsam/bam:\n' + str(sam_seq_ID_list) + '\n')
    else:
        sys.exit(
            '\nERROR: GFF & mapping reference sequence IDs are inconsistent. \n' +
            'SYS.EXIT: Please check and update the sequence IDs in your sequence and gff files so they match up before running again.\ngff:\n' +
            str(gff_seq_ID_list) + '\nsam/bam:\n' + str(sam_seq_ID_list) + '\n' +
            'NOTE: If the sequence ID mismatch is benign e.g. an extra plasmid override by using --gff_force\n')

    print(type(sam_seq_ID_list))
    print(type(gff_seq_ID_list))
    print((sam_seq_ID_list))
    print((gff_seq_ID_list))


def coordinates_to_features_reps(sam_stem, attr_to_columns, condition_label):
    coord_reps_df = pd.read_csv(sam_stem + "_insert_coords.txt", sep='\t', dtype={'ref_name': "str",
                                                                                  'coord': "int64",
                                                                                  'read_grp': "str"})

    coord_counts_reps_df = coord_reps_df.groupby(["ref_name", "coord", "read_grp"]).size().reset_index(
        name=condition_label + '_')

    read_grps = sorted(coord_reps_df.read_grp.unique())

    if len(read_grps) < 3:
        print("Warning: Unable to find >= 3 read groups in fastq data, continuing without replicate insertion counts")
        # returning an empty dataframe
        return pd.DataFrame()

    coord_df_pivot = coord_counts_reps_df.copy(deep=False).pivot_table(index=["ref_name", "coord"],
                                                                       columns=['read_grp'],
                                                                       values=[condition_label + '_'],
                                                                       fill_value=0).reset_index()

    coord_df_pivot.columns = [''.join(col).strip() for col in coord_df_pivot.columns.values]

    old_rep_names = [condition_label + '_' + str(x) for x in read_grps]
    new_rep_names = [condition_label + '_' + "MP" + str(x) for x in range(1, len(read_grps) + 1)]

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
    ## wierd sqlite concatenation + >> || ,  '%' == wildcard double check effect of this
    # this line should allow multi contig files

    mp_coords_join_gff = ps.sqldf(sqlcode, locals())

    # remove first 2 columns ref_name, coord sums the rest according to the feature coordinate groups
    mp_reps_feature_counts = mp_coords_join_gff.drop(mp_coords_join_gff.columns[[0, 1]], axis=1).groupby(
        ['seq_id', 'start', 'end']).agg(["sum"]).reset_index()

    mp_reps_feature_counts.columns = mp_reps_feature_counts.columns.get_level_values(0)
    return mp_reps_feature_counts

    # coord_df_pivot.to_csv("pimms_insert_MP_coordinates_" + condition_label + ".txt", index=False, sep='\t', header=True)
    # attr_to_columns.to_csv("pimms_insert_attr_to_columns_" + condition_label + ".txt", index=False, sep='\t',
    #                        header=True)
    # gff_columns_addback.to_csv("pimms_insert_gff_columns_addback_" + condition_label + ".txt", index=False, sep='\t',
    #                            header=True)


def coordinates_to_features(sam_stem, attr_to_columns, gff_columns_addback, condition_label):
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
    coord_counts_df_pimms2_gff.to_csv("pimms_insert_coordinates_" + condition_label + ".gff", index=False, sep='\t',
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
    ## wierd sqlite concatenation + >> || ,  '%' == wildcard double check effect of this
    # this line should allow multi contig files

    coords_join_gff = ps.sqldf(sqlcode, locals())

    # debugging save of intermediate data
    coords_join_gff.to_csv("pimms_coords_join_gffstuff" + condition_label + ".txt", index=False, sep='\t', header=False)

    # quick dataframe summary
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

    # coords_join_gff.to_csv("pimms_coords_join_gffstuff_test01" + condition_label + ".txt", index=False, sep='\t', header=False)

    # pimms_result_table_test01_notgrouped = coords_join_gff[['seq_id', 'locus_tag', 'gene', 'start', 'end', 'feat_length', 'product'] + gff_extra]
    # pimms_result_table_test01_notgrouped.to_csv("pimms_result_table_test01_notgrouped_" + condition_label + ".txt", index=False, sep='\t', header=True)
    # coords_join_gff_fillna = coords_join_gff.fillna({'a': 0, 'b': 0})
    # coords_join_gff_fillna = \

    ### Important fix group by doesn't work -- any rows with nan values get dropped *yikes very bad!!!!*
    coords_join_gff.fillna('', inplace=True)
    # pimms_result_table_test01 = coords_join_gff.groupby(
    #     ['seq_id', 'locus_tag', 'gene', 'start', 'end', 'feat_length', 'product']).agg(
    #     num_insertions_mapped_per_feat=('counts', 'sum')#,
    #     #num_insert_sites_per_feat=('counts', 'count'),
    #     #first_insert_posn_as_percentile=('posn_as_percentile', 'min'),
    #     #last_insert_posn_as_percentile=('posn_as_percentile', 'max')
    # ).reset_index()
    #
    # pimms_result_table_test01.to_csv("pimms_result_table_test01_grouped_fillna1_" + condition_label + ".txt", index=False, sep='\t', header=True)

    pimms_result_table = coords_join_gff.groupby(
        ['seq_id', 'ID',  # added ID field as unique identifier
         'locus_tag', 'type', 'gene', 'start', 'end', 'feat_length', 'product'] + gff_extra).agg(
        num_insertions_mapped_per_feat=('counts', 'sum'),
        num_insert_sites_per_feat=('counts', 'count'),
        first_insert_posn_as_percentile=('posn_as_percentile', 'min'),
        last_insert_posn_as_percentile=('posn_as_percentile', 'max')
    ).reset_index()

    pimms_result_table.to_csv("pimms_coords_join_prt1_" + condition_label + ".txt", index=False, sep='\t', header=True)

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

    pimms_result_table.to_csv("pimms_coords_join_prt2_" + condition_label + ".txt", index=False, sep='\t', header=False)

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

    NAvalues = {'num_insertions_mapped_per_feat': int(0),
                'num_insert_sites_per_feat': int(0),
                'num_insert_sites_per_feat_per_kb': int(0),
                'first_insert_posn_as_percentile': int(0),
                'last_insert_posn_as_percentile': int(0),
                'NRM_score': int(0),
                'NIM_score': int(0)}
    pimms_result_table_full = pd.merge(gff_columns_addback, pimms_result_table, how='left').fillna(value=NAvalues)

    gff_columns_addback.to_csv("pimms_coords_join_gff_columns_addback_" + condition_label + ".txt", index=False,
                               sep='\t', header=False)
    pimms_result_table_full.to_csv("pimms_coords_join_prtf1_" + condition_label + ".txt", index=False, sep='\t',
                                   header=False)

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


# main function
# def main():

ap = configargparse.ArgumentParser(  # description='PIMMS2 sam/bam processing',
    prog="demo_pimms2_sam_extract",
    add_config_file_help=False,
    config_file_parser_class=configargparse.DefaultConfigFileParser,
    epilog="\n\n*** N.B. This is a development version ***\n \n ",
    description='''description here'''
)
ap.add_argument('--version', action='version', version='%(prog)s 2.0.1 demo')
ap.add_argument('--nano', action='store_true', required=False,
                help='global setting to change processing to nanopore data [False: illumina processing]')
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
samcoords.add_argument("--sam", required=True, nargs=1, metavar='pimms.sam/bam', type=extant_file,
                       help="sam/bam file of mapped IS flanking sequences ")
samcoords.add_argument("--label", required=False, nargs=1, metavar='condition_name', default=[''],
                       help="text tag to add to results file")
samcoords.add_argument("--mismatch", required=False, nargs=1, type=float, metavar='float', default=[None],
                       choices=[Range(0.0, 0.2)],
                       help="fraction of permitted mismatches in mapped read ( 0 <= float < 0.2 [no filter]")
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

# exit()
# find required orog in path or exit with message
# prog_in_path_check('minimap2')

# sort out extra requested gff annotation fields
if parsed_args[0].gff_extra:
    # strip any formatting quotes and turn comma separated string into a list of fields
    gff_extra = parsed_args[0].gff_extra[0].strip("'\"").split(',')
else:
    gff_extra = []

print("extra gff fields: " + str(gff_extra))
gff_file = parsed_args[0].gff[0]
gff_feat_type = ['CDS', 'tRNA', 'rRNA']
# process the gff file to get required fields
gff_columns_addback, attr_to_columns = process_gff(gff_file, gff_feat_type, gff_extra)
# trap multiple return values from function
# gff_file = '/Users/svzaw/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/S.uberis_0140J.gff3'
# gff_file = 'UK15_assembly_fix03.gff'
gff_columns_addback_pseudo, attr_to_columns_pseudo = process_gff(gff_file, ['pseudogene'], [])

# print(gff_columns_addback_pseudo)
# print(attr_to_columns_pseudo)

# exit()
# set up various variables and commandline parameters
use_fraction_mismatch = True
min_depth_cutoff = parsed_args[0].min_depth[0]
# permitted_mismatch = 6
fraction_mismatch = parsed_args[0].mismatch[0]
# sam_file = "0140J_test.IN.PIMMS.RX.rmIS.sub0_min25_max50.ngm_sensitive.90pc.bam"
# sam_file = "UK15_Media_RX_pimms2out_trim60_v_pGh9_UK15.bam"
# sam_file = "UK15_Blood_RX_pimms2out_trim60_v_pGh9_UK15.bam"
sam_file = parsed_args[0].sam[0]
condition_label = parsed_args[0].label[0]

# sam_stem, sam_ext = os.path.splitext(os.path.basename(sam_file))
# sam_stem = sam_stem + '_md' + str(min_depth_cutoff) + '_mm' + str(fraction_mismatch or '')

# check all sequence names match up
seqID_consistancy_check(gff_columns_addback, sam_file)
# process pimms sam/bam  file and produce coordinate / bed files
process_sam(sam_file)

sam_stem = modify_sam_stem(sam_file)

# allocate insertions to features and create results merged with GFF
# possibly poor coding to merge with gff here
pimms_result_table_full = coordinates_to_features(sam_stem, attr_to_columns, gff_columns_addback, condition_label)

if not parsed_args[0].noreps:
    mp_reps_feature_counts = coordinates_to_features_reps(sam_stem, attr_to_columns, condition_label)
    if not mp_reps_feature_counts.empty:
        merged_with_reps = pimms_result_table_full.merge(mp_reps_feature_counts, on=["seq_id", "start", "end"],
                                                         how='inner')
        pimms_result_table_full = merged_with_reps

if not gff_columns_addback_pseudo.empty:
    tag_psueudogenes = gff_columns_addback_pseudo['locus_tag']
    pimms_result_table_full.loc[pimms_result_table_full.locus_tag.isin(tag_psueudogenes), "type"] = \
        pimms_result_table_full['type'] + '_pseudo'

# write results as text/excel
pimms_result_table_full.to_csv(sam_stem + "_countinfo_tab.txt", index=False, sep='\t')
pimms_result_table_full.to_csv(sam_stem + "_countinfo.csv", index=False, sep=',')
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


#    if __name__ == '__main__':
#       main()
