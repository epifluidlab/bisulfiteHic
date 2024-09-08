import numpy as np
from scipy.stats import pearsonr, spearmanr
import argparse
import pyBigWig
import multiprocessing
from functools import partial

def process_line(bw_tracks, line):
    line = line.rstrip('\n')
    chrom1, start1, end1, chrom2, start2, end2, name, _ = line.split("\t")
    methy_locus1 = []
    methy_locus2 = []
    for bw in bw_tracks:
        stat1 = bw.stats(chrom1, int(start1), int(end1), exact=True)[0]
        stat2 = bw.stats(chrom2, int(start2), int(end2), exact=True)[0]
        if stat1 is not None and stat2 is not None:
            methy_locus1.append(stat1)
            methy_locus2.append(stat2)
    if len(methy_locus1) > 0:
        methy_locus1 = np.array(methy_locus1, dtype=float)
        methy_locus2 = np.array(methy_locus2, dtype=float)
        cor_val, p_val = pearsonr(methy_locus1, methy_locus2)
        spearman_cor_val, spearman_p_val = spearmanr(methy_locus1, methy_locus2)
        num_point = len(methy_locus1)
        return f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{name}\t{cor_val:.3f}\t{p_val:.3f}\t{spearman_cor_val:.3f}\t{spearman_p_val:.3f}\t{num_point}\n"
    else:
        return None

def cal_concord_sc_at_regions(methy_file_list, bedpe, out_summary, threads=1):
    bw_tracks = [pyBigWig.open(line.strip()) for line in open(methy_file_list, 'r')]

    # Open output file
    with open(out_summary, 'w') as out:
        # Read all lines from the bedpe file into memory
        lines = open(bedpe, 'r').readlines()

        # Create a pool of workers
        with multiprocessing.Pool(processes=threads) as pool:
            # Use functools.partial to create a function with some arguments pre-filled
            func = partial(process_line, bw_tracks)

            # Map func over the lines using pool.map
            results = pool.map(func, lines)

            # Write results that are not None to the output file
            for result in results:
                if result:
                    out.write(result)

    # Close all BigWig tracks
    for bw in bw_tracks:
        bw.close()

##output loop bedpe + pearson cor, p value, spearman cor, p value, num point in anchors
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--methy_file_list') #the big wig file list for methylation level (GCH or HCG or CG)
    p.add_argument('--bedpe') #a bedpe file that need to summarize the correlation level
    p.add_argument('--out_summary') ##output the summary of pearson correlation and spearman correlation, each row represent a correlation across cells at each paired region
    p.add_argument('--threads', type=int, default=1)
    args = p.parse_args()
    cal_concord_sc_at_regions(**vars(args))