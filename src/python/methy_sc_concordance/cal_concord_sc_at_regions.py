import numpy as np
from scipy.stats import pearsonr, spearmanr
import argparse
import pyBigWig
import multiprocessing
from functools import partial
from tqdm import tqdm

def process_line(methy_file_list, line):
    bw_tracks = [pyBigWig.open(filename) for filename in methy_file_list]
    elements = line.split("\t")
    # Ensure there are at least 6 elements to unpack
    if line.startswith("#") or len(elements) < 6:
        return None  # Or handle error more specifically if needed
    chrom1, start1, end1, chrom2, start2, end2 = elements[:6]
    try:
        start1 = int(start1)
        end1 = int(end1)
        start2 = int(start2)
        end2 = int(end2)
    except ValueError:
        return f"Error converting positions to integers in line: {line}\n"

    methy_locus1 = []
    methy_locus2 = []
    for bw in bw_tracks:
        chrom1_len = bw.chroms(chrom1)
        chrom2_len = bw.chroms(chrom2)
        stat1 = None
        stat2 = None
        if chrom1_len and chrom2_len and int(end1) < chrom1_len and int(end2) < chrom2_len:
            stat1 = bw.stats(chrom1, int(start1), int(end1), exact=True)[0]
            stat2 = bw.stats(chrom2, int(start2), int(end2), exact=True)[0]
        else:
            if chrom1_len and int(end1) >= chrom1_len:
                stat1 = bw.stats(chrom1, int(start1), chrom1_len, exact=True)[0]
            if chrom2_len and int(end2) >= chrom2_len:
                stat2 = bw.stats(chrom2, int(start2), chrom2_len, exact=True)[0]
        if stat1 is not None and stat2 is not None:
                methy_locus1.append(stat1)
                methy_locus2.append(stat2)
        bw.close()
    methy_locus1 = np.array(methy_locus1, dtype=float)
    methy_locus2 = np.array(methy_locus2, dtype=float)
    #methy_locus1 = methy_locus1[~np.isnan(methy_locus1) & ~np.isinf(methy_locus1)]
    #methy_locus2 = methy_locus2[~np.isnan(methy_locus2) & ~np.isinf(methy_locus2)]
    if methy_locus1.size > 10:
        cor_val, p_val = pearsonr(methy_locus1, methy_locus2)
        spearman_cor_val, spearman_p_val = spearmanr(methy_locus1, methy_locus2)
        num_point = len(methy_locus1)
        name = f"{chrom1}:{start1}:{end1}:{chrom2}:{start2}:{end2}"
        return f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{name}\t{cor_val:.3f}\t{p_val:.3f}\t{spearman_cor_val:.3f}\t{spearman_p_val:.3f}\t{num_point}\n"
    return None

def cal_concord_sc_at_regions(methy_file_list, bedpe, out_summary, threads=1):
    # Read filenames into a list
    with open(methy_file_list, 'r') as f:
        filenames = [line.strip() for line in f.readlines()]

    # Open the output file
    with open(out_summary, 'w') as out:
        # Read all lines from the bedpe file into memory
        lines = open(bedpe, 'r').readlines()

        # Create a pool of workers
        with multiprocessing.Pool(processes=threads) as pool:
            func = partial(process_line, filenames)
            results = pool.imap(func, lines)

            # Use tqdm for progress update
            for result in tqdm(results, total=len(lines), desc="Processing lines"):
                if result:
                    out.write(result)


##output loop bedpe + pearson cor, p value, spearman cor, p value, num point in anchors
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--methy_file_list') #the big wig file list for methylation level (GCH or HCG or CG)
    p.add_argument('--bedpe') #a bedpe file that need to summarize the correlation level
    p.add_argument('--out_summary') ##output the summary of pearson correlation and spearman correlation, each row represent a correlation across cells at each paired region
    p.add_argument('--threads', type=int, default=1)
    args = p.parse_args()
    cal_concord_sc_at_regions(**vars(args))