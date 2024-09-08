import numpy as np
from scipy.stats import pearsonr, spearmanr
import argparse
import pyBigWig

def cal_concord_sc_at_regions(methy_file_list, bedpe, out_summary):
    bw_tracks = []
    with open(methy_file_list, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            bw_tracks.append(pyBigWig.open(line))

    with open(out_summary, 'w') as out:
        with open(bedpe, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                chrom1, start1, end1, chrom2, start2, end2, name, _ = line.split("\t")
                methy_locus1 = []
                methy_locus2 = []
                for i in range(0, len(bw_tracks)):
                    bw = bw_tracks[i]
                    stat1 = bw(chrom1, int(start1), int(end1), exact=True)[0]
                    stat2 = bw(chrom2, int(start2), int(end2), exact=True)[0]
                    if stat1 is not None and stat2 is not None:
                        methy_locus1.append(stat1)
                        methy_locus2.append(stat2)
                if len(methy_locus1) > 0:
                    methy_locus1 = np.array(methy_locus1,type="float")
                    methy_locus2 = np.array(methy_locus2,type="float")
                    cor_val, p_val = pearsonr(methy_locus1, methy_locus2)
                    spearman_cor_val, spearman_p_val = spearmanr(methy_locus1, methy_locus2)
                    num_point = len(methy_locus1)
                    out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%d\n" % (chrom1, start1, end1, chrom2, start2, end2, name, cor_val, p_val, spearman_cor_val, spearman_p_val, num_point))

    for bw in bw_tracks:
        bw.close()

##output loop bedpe + pearson cor, p value, spearman cor, p value, num point in anchors
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--methy_file_list') #the big wig file list for methylation level (GCH or HCG or CG)
    p.add_argument('--bedpe') #a bedpe file that need to summarize the correlation level
    p.add_argument('--out_summary') ##output the summary of pearson correlation and spearman correlation, each row represent a correlation across cells at each paired region
    args = p.parse_args()
    cal_concord_sc_at_regions(**vars(args))