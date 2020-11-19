import numpy as np
import argparse
from bx.intervals.intersection import Intersecter, Interval


def get_gene_density(file = '/home/dnaase/startup/data/annotations/gene/gencode.v30.b37.gene_density.bedgraph',
                            chrom_len_dict = None, chrom=None,
                            w_len=2500000):
    genes_dict = dict()
    with open(file,'r') as f:
        for line in f:
            line=line.rstrip('\n')
            chrom_i,start,end,_ = line.split("\t")
            if chrom_i in genes_dict:
                genes = genes_dict[chrom_i]
            else:
                genes = Intersecter()
            genes.add_interval(Interval(int(start), int(end), chrom=chrom_i))
            genes_dict[chrom_i] = genes
    genes = genes_dict[chrom]
    gene_density = []
    for w_start in range(0,chrom_len_dict[chrom],w_len):
        w_end = min(chrom_len_dict[chrom], w_start + w_len)
        gene_intersect = genes.find(w_start, w_end)
        gene_density.append(len(gene_intersect)*1000000/w_len)
    return np.array(gene_density,dtype='float')


def eigen_adjust_direction(in_eig_file, out_file, chrom_size, w_len, gene_density):
    chrom_len_dict = dict()
    with open(chrom_size, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            chrom, size = line.split("\t")
            chrom_len_dict[chrom] = int(size)
    gene_density = get_gene_density(file=gene_density, chrom_len_dict = chrom_len_dict, w_len=w_len)
    map_1D_mat = []

    with open(in_eig_file,'r') as f:
        for line in f:
            eig=line.rstrip('\n')
            map_1D_mat.append(eig)
    map_1D_mat = np.array(map_1D_mat, dtype='float')
    no_nan_index = [i for i, s in enumerate(map_1D_mat) if (not np.isnan(s))]
    gene_density = gene_density[no_nan_index]
    map_1D_mat = map_1D_mat[no_nan_index]
    pc1_pos_index = [i for i, s in enumerate(map_1D_mat) if (s >= 0)]
    pc1_neg_index = [i for i, s in enumerate(map_1D_mat) if (s < 0)]
    if np.nanmean(gene_density[pc1_pos_index]) < np.nanmean(gene_density[pc1_neg_index]):
        map_1D_mat = 0 - map_1D_mat
    w_start = 0
    with open(out_file, 'w') as out:
        for i in range(len(map_1D_mat)):
            w_end = min(chrom_len_dict[chrom], w_start + w_len)
            out.write(chrom + '\t' + str(w_start) + '\t' + str(w_end) + '\t' + str(map_1D_mat[i]) + '\n')
            w_start = w_start + w_len


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--in_eig_file') # the eigenvector file from juicer tools, just the eigen value for that chromosome.
    p.add_argument('--out_file') # chr, start, end, adjusted eigen value.
    p.add_argument('--chr') ##specify which chr
    p.add_argument('--chrom_size', default='/home/dnaase/startup/data/genomes/hg19/human_g1k_v37.common_chrom.size')
    p.add_argument('--w_len', type=int, default=2500000)
    p.add_argument('--gene_density',default = '/home/dnaase/startup/data/annotations/gene/gencode.v30.b37.gene_density.bedgraph')
    args = p.parse_args()
    eigen_adjust_direction(**vars(args))
