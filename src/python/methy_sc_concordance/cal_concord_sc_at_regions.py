import sys, glob, gzip, subprocess, re, os
import numpy as np
import collections
import argparse
import pyBigWig

def cal_concord_sc_across_cells_per_region(chrom1, start1, end1, chrom2, start2, end2, bw_methy, bw_cov, min_point=1):



def cal_concord_sc_at_regions(methy_file_list, cov_file_list, bedpe, out_summary, min_point=1):
    bw_tracks = []
    with open(methy_file_list, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            bw_tracks.append(pyBigWig.open(line))

    bw_cov_tracks = []
    with open(cov_file_list, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            bw_cov_tracks.append(pyBigWig.open(line))

    with open(out_summary, 'w') as out:
        with open(bedpe, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                chrom1, start1, end1, chrom2, start2, end2, name, _ = line.split("\t")



    for bw in bw_tracks:
        bw.close()
    for bw in bw_cov_tracks:
        bw.close()

    meta_track = pyBigWig.open(file)
    meta_row = []
    for w_start in gene_positions:
        w_end = min(chrom_len_dict[chrom], w_start + w_len)
        if coverage:
            meta_row.append(meta_track.stats(chrom, w_start, w_end, type='coverage', exact=True))
        else:
            meta_row.append(meta_track.stats(chrom, w_start, w_end, exact=True))
    meta_row = np.array(meta_row, dtype='float').reshape(-1)
    meta_track.close()


# due to the error of import pyBigWig and bx-python in mac OSX
def wrap_bigwig_tool(chrom, start, end):
    result = subprocess.run(['ls', (chrom + " " + str(start) + " " + str(end) + " " + str(1))], stdout=subprocess.PIPE)
    result = result.stdout.decode('utf-8')
    if result:
        return result
    else:
        return np.nan



##output loop bedpe + pearson cor, p value, spearman cor, p value, num point in anchor 1
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--methy_file_list') #the big wig file list for methylation level (GCH or HCG)
    p.add_argument('--cov_file_list') #the big wig file list for total coverage level (GCH or HCG), should be the same order as in --methy_file_list
    p.add_argument('--bedpe') #a bedpe file that need to summarize the correlation level
    p.add_argument('--out_summary') ##output the summary of pearson correlation and spearman correlation, each row represent a correlation across cells at each paired region
    p.add_argument('--min_point', type=int, default=1) ## minimum number of points at each paired region.
    args = p.parse_args()
    cal_concord_sc_at_regions(**vars(args))