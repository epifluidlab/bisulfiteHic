import argparse, gzip
from bx.intervals.intersection import Intersecter, Interval
from scipy import stats
import numpy as np
import pysam


def get_alpha_beta(mu, var):
    if var == 0 or mu == 0:
        alpha = np.nan
        beta = np.nan
    else:
        alpha = ((1 - mu) / var - 1 / mu) * mu * mu
        beta = alpha * (1 / mu - 1)
    return (alpha, beta)


def betabinom_test_mgch_per_chrom(tbx, out, chrom, size, min_gch, min_point, window1, window2):
    intersecter = Intersecter()
    num = 0
    c_global = 0
    m_global = 0
    global_methy = []
    for row in tbx.fetch(chrom, 0, size):
        chrom, start, end, _, score, strand, m, c = row.split("\t")
        interv_tmp = Interval(int(start), int(end), chrom=chrom, strand=strand,
                              value={'score': int(score), 'methylation': float(m), 'coverage': int(c)})
        intersecter.add_interval(interv_tmp)
        c_global = c_global + int(c)
        m_global = m_global + int(float(m) * int(c) / 100.0)
        global_methy.append(float(m) / 100.0)
        num = num + 1
    print('with size of GCH:' + str(num))
    if num < 10000:
        print('too small number of GCH in this chrom, skip')
        return
    global_methy = np.array(global_methy, dtype='float32')
    global_methy_mu = np.nanmean(global_methy)
    global_methy_var = np.nanvar(global_methy)
    global_methy_alpha, global_methy_beta = get_alpha_beta(global_methy_mu, global_methy_var)
    print('global alpha:' + str(global_methy_alpha))
    print('global beta:' + str(global_methy_beta))

    del global_methy

    for line in tbx.fetch(chrom, 0, size):
        chrom, start, end, _, score, strand, m, c = line.split("\t")
        methy_count = int(float(m) * int(c) / 100.0)
        total_count = int(c)
        methy = float(m) / 100
        num_gch_1 = 0
        c_window_1 = 0
        m_window_1 = 0
        meth_list1 = []
        w1_start_1 = max(int(start) - window1, 0)
        w1_end_1 = int(start)
        w1_start_2 = int(end)
        w1_end_2 = min(int(end) + window1, size)
        for interval in intersecter.find(w1_start_1, w1_end_1):
            c_int = interval.value['coverage']
            c_window_1 = c_window_1 + c_int
            m_window_1 = m_window_1 + int(interval.value['methylation'] * c_int / 100.0)
            meth_list1.append(interval.value['methylation'] / 100.0)
            num_gch_1 = num_gch_1 + 1
        for interval in intersecter.find(w1_start_2, w1_end_2):
            c_int = interval.value['coverage']
            c_window_1 = c_window_1 + c_int
            m_window_1 = m_window_1 + int(interval.value['methylation'] * c_int / 100.0)
            meth_list1.append(interval.value['methylation'] / 100.0)
            num_gch_1 = num_gch_1 + 1
        meth_list1 = np.array(meth_list1, dtype='float32')
        local_methy1_mu = np.nanmean(meth_list1)
        local_methy1_var = np.nanvar(meth_list1)
        local_methy1_alpha, local_methy1_beta = get_alpha_beta(local_methy1_mu, local_methy1_var)
        if c_window_1 > 0:
            local_fold1 = methy - (m_window_1 / c_window_1)
        else:
            local_fold1 = methy - m_global/c_global

        num_gch_2 = 0
        c_window_2 = 0
        m_window_2 = 0
        meth_list2 = []
        w2_start_1 = max(int(start) - window2, 0)
        w2_end_1 = int(start)
        w2_start_2 = int(end)
        w2_end_2 = min(int(end) + window2, size)
        for interval in intersecter.find(w2_start_1, w2_end_1):
            c_int = interval.value['coverage']
            c_window_2 = c_window_2 + c_int
            m_window_2 = m_window_2 + int(interval.value['methylation'] * c_int / 100.0)
            meth_list2.append(interval.value['methylation'] / 100.0)
            num_gch_2 = num_gch_2 + 1
        for interval in intersecter.find(w2_start_2, w2_end_2):
            c_int = interval.value['coverage']
            c_window_2 = c_window_2 + c_int
            m_window_2 = m_window_2 + int(interval.value['methylation'] * c_int / 100.0)
            meth_list2.append(interval.value['methylation'] / 100.0)
            num_gch_2 = num_gch_2 + 1
        meth_list2 = np.array(meth_list2, dtype='float32')
        local_methy2_mu = np.nanmean(meth_list2)
        local_methy2_var = np.nanvar(meth_list2)
        local_methy2_alpha, local_methy2_beta = get_alpha_beta(local_methy2_mu, local_methy2_var)
        if c_window_2 > 0:
            local_fold2 = methy - (m_window_2 / c_window_2)
        else:
            local_fold2 = methy - (m_global / c_global)
        global_fold = methy - (m_global / c_global)

        if num_gch_1 >= min_gch and c_window_1 >= min_point and (not np.isnan(local_methy1_alpha)) and (not np.isnan(local_methy1_beta)):
            window1_p = m_window_1 / c_window_1
            beta_binom_p_w1 = 0 - np.log10(stats.betabinom.cdf(methy_count, total_count, local_methy1_alpha,
                                                         local_methy1_beta))
        else:
            window1_p = m_global / c_global
            beta_binom_p_w1 = 0 - np.log10(stats.betabinom.cdf(methy_count, total_count, global_methy_alpha, global_methy_beta))
        binomial_p_w1 = 0 - np.log10(stats.binom_test(methy_count, n=total_count, p=window1_p, alternative='less'))
        if num_gch_2 >= min_gch and c_window_2 >= min_point and (not np.isnan(local_methy2_alpha)) and (not np.isnan(local_methy2_beta)):
            window2_p = m_window_2 / c_window_2
            beta_binom_p_w2 = 0 - np.log10(stats.betabinom.cdf(methy_count, total_count, local_methy2_alpha, local_methy2_beta))
        else:
            window2_p = m_global / c_global
            beta_binom_p_w2 = 0 - np.log10(stats.betabinom.cdf(methy_count, total_count, global_methy_alpha, global_methy_beta))
        binomial_p_w2 = 0 - np.log(stats.binom_test(methy_count, n=total_count, p=window2_p, alternative='less'))
        binomial_p_global = 0 - np.log10(
            stats.binom_test(methy_count, n=total_count, p=m_global / c_global, alternative='less'))
        binomial_p_final = np.nanmax((binomial_p_w1, binomial_p_w2, binomial_p_global))
        beta_binom_p_global = 0 - np.log10(stats.betabinom.cdf(methy_count, total_count, global_methy_alpha, global_methy_beta))
        beta_binom_p_final = np.nanmax((beta_binom_p_w1, beta_binom_p_w2, beta_binom_p_global))
        out.write(line + '\t' + str(round(local_fold1, 4)) + '\t' + str(round(local_fold2, 4)) + '\t' + str(
            round(global_fold, 4)) + '\t' + str(round(binomial_p_w1, 4)) + '\t' + str(
            round(binomial_p_w2, 4)) + '\t' + str(round(binomial_p_global, 4)) + '\t' + str(
            binomial_p_final) + '\t' + str(round(beta_binom_p_w1, 4)) + '\t' + str(
            round(beta_binom_p_w2, 4)) + '\t' + str(round(beta_binom_p_global, 4)) + '\t' + str(
            beta_binom_p_final) + '\n')


def betabinom_test_mgch(in_bed, out_bed, chrom_size, min_gch, min_point, window1, window2):
    chrom_dict = dict()
    with open(chrom_size, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            splitline = line.split("\t")
            chrom = splitline[0]
            size = int(splitline[1])
            chrom_dict[chrom] = size
    tbx = pysam.TabixFile(in_bed)
    if out_bed.endswith('.gz'):
        out = gzip.open(out_bed, 'wt')
    else:
        out = open(out_bed, 'w')
    for chrom in chrom_dict.keys():
        if chrom in tbx.contigs:
            print('chrom:' + chrom)
            size = chrom_dict[chrom]
            betabinom_test_mgch_per_chrom(tbx, out, chrom, size, min_gch, min_point, window1, window2)
    out.close()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--in_bed')  ##bgzip and indexed 6plus2 bed file, standard bed + methy (%, 0-100 scale) + coverage
    p.add_argument(
        '--out_bed')  ## standard input bed, then log2 fold change (local bg 1, local bg 2, global bg, min), -log10p (local bg 1, local bg 2, global bg, max)
    p.add_argument('--chrom_size')  ## to get accessibility background in each chr; recommend only autosome + chrX
    p.add_argument('--window1', type=int, default=2500)  ##left 2500bp and right 2500bp
    p.add_argument('--window2', type=int, default=5000)  ##left 5000bp and right 5000bp
    p.add_argument('--min_gch', type=int, default=10)  ##requires at least 10 GCH
    p.add_argument('--min_point', type=int,
                   default=100)  ##requires at least 100 data points
    args = p.parse_args()
    betabinom_test_mgch(**vars(
        args))  ##get beta-binomial test, -log10 pvalue, for regions greater methy, "greater", regions smaller methy, use "smaller"hypothesis testing. # maybe also get FDR value
