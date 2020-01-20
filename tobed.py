import os
import datetime
import multiprocessing
import pysam
import pandas as pd
import gffpandas.gffpandas as gffpd
import pandasql as ps
import pysam


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
def main():
    insamfile = "0140J_test.IN.PIMMS.RX.rmIS.sub0_min25_max50.ngm_sensitive.90pc.bam"
    annotation = gffpd.read_gff3('/Users/svzaw/Data/PIMMS_redo/PIMMS2_stuff/PIMMS_V1/S.uberis_0140J.gff3')
    attr_to_columns = annotation.attributes_to_columns()
    attr_to_columns = attr_to_columns[attr_to_columns['type'] == 'gene']
    attr_to_columns = attr_to_columns.assign(feat_length=(attr_to_columns.end - attr_to_columns.start + 1)).drop(
        columns=['attributes'])

    samfile = pysam.AlignmentFile(insamfile, "rb")
    open("test.bed", 'w').close()
    f = open("test.bed", "a")

    STRAND = ["+", "-"]
    for read in samfile.fetch():
        STR = STRAND[int(read.is_reverse)]
        BED = [read.reference_name, read.pos, read.reference_end, ".", read.mapq, STR, '# ' + read.query_name]
        f.write('\t'.join([str(i) for i in BED]))
        f.write('\n')
    f.close()

    open("test.insert_coords.txt", 'w').close()
    f2 = open("test.insert_coords.txt", "a")
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

    coord_df = pd.read_csv("test.insert_coords.txt", sep='\t', dtype={'coord': "int64"})
    coord_counts_df = coord_df.groupby(['ref_name', 'coord']).size().reset_index(name='counts')
    coord_counts_gt1_df = coord_counts_df[coord_counts_df['counts'] > 1]

    sqlcode = '''
    select coord_counts_df.*
    ,attr_to_columns.*
    from coord_counts_df
    inner join attr_to_columns 
    on coord_counts_df.coord between attr_to_columns.start and attr_to_columns.end
    '''

    coords_join_gff = ps.sqldf(sqlcode, locals())

    ## add position as percentile (needs manual confirmation)
    coords_join_gff = coords_join_gff.assign(
        posn_as_percentile=(((coords_join_gff.coord - coords_join_gff.start) + 1) / (
                coords_join_gff.feat_length / 100)).where(
            coords_join_gff.strand == '+', ((coords_join_gff.end - coords_join_gff.coord) + 1) / (
                    coords_join_gff.feat_length / 100))).round(
        {"posn_as_percentile": 1})  # =(attr_to_columns.end - attr_to_columns.start + 1))
    #  sprintf("%.1f", ((($in-$in_start)+1) / ($in_length/100)));
    #  sprintf("%.1f", ((($in_stop-$ in) + 1) / ($in_length / 100)));

    #    ,attr_to_columns.start
    #    ,attr_to_columns.end
    # coord_counts_df.coord
    #     ,coord_counts_df.ref_name
    #     ,coord_counts_df.counts
    #     ,* from attr_to_columns
    # sqlcode = '''
    # select df_1.timestamp
    # ,df_1.A
    # ,df_1.B
    # ,df_2.event
    # from df_1
    # inner join df_2
    # on d1.timestamp between df_2.start and df2.end
    # '''
    #
    # newdf = ps.sqldf(sqlcode,locals())

    if __name__ == '__main__':
        main()
