import os
import datetime
import multiprocessing
import pysam
import pandas as pd
import gffpandas.gffpandas as gffpd
import pandasql as ps
import pysam
import urllib as ul

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


# main function
# def main():

min_depth_cutoff = 1
# insamfile = "0140J_test.IN.PIMMS.RX.rmIS.sub0_min25_max50.ngm_sensitive.90pc.bam"
insamfile = "pGh9_UK15_UK15_Media_Input_pimmsout_trim50_nodecon_mm2.bam"
# annotation = gffpd.read_gff3('/Users/svzaw/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/S.uberis_0140J.gff3')
annotation = gffpd.read_gff3('UK15_genome_fix02.gff')
attr_to_columns = annotation.attributes_to_columns()
attr_to_columns = attr_to_columns[attr_to_columns['type'] == 'CDS']
attr_to_columns = attr_to_columns.assign(feat_length=(attr_to_columns.end - attr_to_columns.start + 1)).dropna(axis=1,
                                                                                                               how='all').drop(
    columns=['attributes'])

# remove RFC 3986 % encoding from product (gff3 attribute)
attr_to_columns = attr_to_columns.assign(product_nopc=attr_to_columns['product'].apply(ul.parse.unquote)).drop(
    columns=['product']).rename(columns={'product_nopc': 'product'})

# attr_to_columns = attr_to_columns.dropna(axis=1, how='all')
samfile = pysam.AlignmentFile(insamfile, "rb")
open("test2.bed", 'w').close()
f = open("test2.bed", "a")

STRAND = ["+", "-"]
for read in samfile.fetch():
    STR = STRAND[int(read.is_reverse)]
    BED = [read.reference_name, read.pos, read.reference_end, ".", read.mapq, STR, '# ' + read.query_name]
    f.write('\t'.join([str(i) for i in BED]))
    f.write('\n')
f.close()

open("test2.insert_coords.txt", 'w').close()
f2 = open("test2.insert_coords.txt", "a")
f2.write('\t'.join([str(i) for i in ['ref_name', 'coord', 'strand', 'read_name']]))
f2.write('\n')
STRAND = ["+", "-"]
for read in samfile.fetch():
    STR = STRAND[int(read.is_reverse)]
    if STR == '+':
        COORDS = [read.reference_name, (read.pos + 4), STR, '# ' + read.query_name]
        f2.write('\t'.join([str(i) for i in COORDS]))
        f2.write('\n')
    if STR == '-':
        COORDS = [read.reference_name, (read.reference_end - 4), STR, '# ' + read.query_name]
        f2.write('\t'.join([str(i) for i in COORDS]))
        f2.write('\n')
f2.close()

coord_df = pd.read_csv("test2.insert_coords.txt", sep='\t', dtype={'coord': "int64"})
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
    from coord_counts_df
    inner join attr_to_columns 
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
    ['seq_id', 'locus_tag', 'gene', 'start', 'end', 'feat_length', 'product']).agg(
    num_reads_mapped_per_feat=('counts', 'sum'),
    num_insert_sites_per_feat=('counts', 'count'),
    posn_as_percentile_first=('posn_as_percentile', 'min'),
    posn_as_percentile_last=('posn_as_percentile', 'max')
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
                                         'product',
                                         'num_reads_mapped_per_feat',
                                         'num_insert_sites_per_feat',
                                         'num_insert_sites_per_feat_per_kb',
                                         'posn_as_percentile_first',
                                         'posn_as_percentile_last',
                                         'NRM_score',
                                         'NIM_score']]

print(list(pimms_result_table.columns.values))

pimms_result_table.to_csv("pimms2_result_table_test2_tab.csv", index=False, sep='\t')

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
