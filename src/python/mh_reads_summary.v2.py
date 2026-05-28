import pysam
import argparse


def mh_reads_summary_v2(in_cram, out_summary):
    total_fragments = 0
    uniq_mapped = 0 ##both end unqiue mapped
    uniq_mapped_mapq30 = 0
    uniq_mapped_mapq30_to_lambda = 0
    uniq_mapped_mapq30_to_target_species = 0
    uniq_mapped_mapq30_nopcr = 0
    uniq_mapped_mapq30_nopcr_to_lambda = 0
    uniq_mapped_mapq30_nopcr_to_target_species = 0
    uniq_mapped_mapq30_nopcr_trans = 0
    uniq_mapped_mapq30_nopcr_cis = 0
    uniq_mapped_mapq30_nopcr_cis_1kb = 0
    uniq_mapped_mapq30_nopcr_cis_10kb = 0
    uniq_mapped_mapq30_nopcr_cis_20kb = 0
    pre_read=None
    with pysam.AlignmentFile(in_cram, "rc") as samfile:
        for read in samfile.fetch(until_eof=True):
            if pre_read is None:
                pre_read = read
                continue
            elif pre_read.query_name != read.query_name:
                pre_read = read
                continue
            elif pre_read.query_name == read.query_name:
                total_fragments = total_fragments + 1
                if read.is_paired and (not read.is_unmapped) and (not pre_read.is_unmapped) and (not read.is_secondary) and (not pre_read.is_secondary):
                    uniq_mapped = uniq_mapped + 1
                    if read.mapping_quality >= 30 and pre_read.mapping_quality >= 30:
                        uniq_mapped_mapq30 = uniq_mapped_mapq30 + 1
                        if read.reference_name == "lambda":
                            uniq_mapped_mapq30_to_lambda = uniq_mapped_mapq30_to_lambda + 1
                        else:
                            uniq_mapped_mapq30_to_target_species = uniq_mapped_mapq30_to_target_species + 1
                        if not read.is_duplicate:
                            uniq_mapped_mapq30_nopcr = uniq_mapped_mapq30_nopcr + 1
                            if read.reference_name == "lambda":
                                uniq_mapped_mapq30_nopcr_to_lambda = uniq_mapped_mapq30_nopcr_to_lambda + 1
                            else:
                                uniq_mapped_mapq30_nopcr_to_target_species = uniq_mapped_mapq30_nopcr_to_target_species + 1
                                if read.reference_name != read.next_reference_name:
                                    uniq_mapped_mapq30_nopcr_trans = uniq_mapped_mapq30_nopcr_trans + 1
                                else:
                                    uniq_mapped_mapq30_nopcr_cis =  uniq_mapped_mapq30_nopcr_cis + 1
                                    # Position-based distance (NOT BAM TLEN).
                                    # bhmem sets TLEN=0 for chimeric Hi-C ligations,
                                    # so abs(read.template_length) under-counts long-range cis.
                                    dist = abs(read.reference_start - pre_read.reference_start)
                                    if dist >= 1000:
                                        uniq_mapped_mapq30_nopcr_cis_1kb = uniq_mapped_mapq30_nopcr_cis_1kb + 1
                                        if dist >= 10000:
                                            uniq_mapped_mapq30_nopcr_cis_10kb = uniq_mapped_mapq30_nopcr_cis_10kb + 1
                                            if dist >= 20000:
                                                uniq_mapped_mapq30_nopcr_cis_20kb = uniq_mapped_mapq30_nopcr_cis_20kb + 1
            else:
                print("Error!")
                break

    # Safe percentage: returns 0.0 if denominator is 0 (avoids ZeroDivisionError
    # for cells with empty or very poor BAMs that have no MAPQ30/no-PCR reads).
    def pct(num, den):
        return (num * 100.0 / den) if den > 0 else 0.0

    with open(out_summary, 'w') as out:
        line = "TotalFragments:\t" + str(total_fragments) + "\t100%" \
               + "\nUniqMapped:\t" + str(uniq_mapped) + "\t" + str(pct(uniq_mapped, total_fragments)) + "%" \
               + "\nUniqMappedMapQ30:\t" + str(uniq_mapped_mapq30) + "\t" + str(pct(uniq_mapped_mapq30, total_fragments)) + "%" \
               + "\nUniqMappedMapQ30ToLambda:\t" + str(uniq_mapped_mapq30_to_lambda) + "\t" + str(
            pct(uniq_mapped_mapq30_to_lambda, uniq_mapped_mapq30)) + "%" \
               + "\nUniqMappedMapQ30ToTargetSpecies:\t" + str(
            uniq_mapped_mapq30_to_target_species) + "\t" + str(
            pct(uniq_mapped_mapq30_to_target_species, uniq_mapped_mapq30)) + "%" \
               + "\nUniqMappedMapQ30NoPcr:\t" + str(uniq_mapped_mapq30_nopcr) + "\t" + str(
            pct(uniq_mapped_mapq30_nopcr, uniq_mapped_mapq30)) + "%" \
               + "\nUniqMappedMapQ30NoPcrToLambda:\t" + str(uniq_mapped_mapq30_nopcr_to_lambda) + "\t" + str(
            pct(uniq_mapped_mapq30_nopcr_to_lambda, uniq_mapped_mapq30)) + "%" \
               + "\nUniqMappedMapQ30NoPcrToTargetSpecies:\t" + str(uniq_mapped_mapq30_nopcr_to_target_species) + "\t" + str(
            pct(uniq_mapped_mapq30_nopcr_to_target_species, uniq_mapped_mapq30)) + "%" \
               + "\nUniqMappedMapQ30NoPcrTrans:\t" + str(uniq_mapped_mapq30_nopcr_trans) + "\t" + str(
            pct(uniq_mapped_mapq30_nopcr_trans, uniq_mapped_mapq30_nopcr_to_target_species)) + "%" \
               + "\nUniqMappedMapQ30NoPcrCis:\t" + str(uniq_mapped_mapq30_nopcr_cis) + "\t" + str(
            pct(uniq_mapped_mapq30_nopcr_cis, uniq_mapped_mapq30_nopcr_to_target_species)) + "%" \
               + "\nUniqMappedMapQ30NoPcrCisMore1kb:\t" + str(uniq_mapped_mapq30_nopcr_cis_1kb) + "\t" + str(
            pct(uniq_mapped_mapq30_nopcr_cis_1kb, uniq_mapped_mapq30_nopcr_to_target_species)) + "%" \
               + "\nUniqMappedMapQ30NoPcrCisMore10kb:\t" + str(uniq_mapped_mapq30_nopcr_cis_10kb) + "\t" + str(
            pct(uniq_mapped_mapq30_nopcr_cis_10kb, uniq_mapped_mapq30_nopcr_to_target_species)) + "%" \
               + "\nUniqMappedMapQ30NoPcrCisMore20kb:\t" + str(uniq_mapped_mapq30_nopcr_cis_20kb) + "\t" + str(
            pct(uniq_mapped_mapq30_nopcr_cis_20kb, uniq_mapped_mapq30_nopcr_to_target_species)) + "%\n"
        out.write(line)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--in_cram') ##required sorted by read_name. R1 and R2's read name should be the same.  not with /1 or /2 suffix.
    p.add_argument('--out_summary') ##only check pair-end reads. only use unique mapped reads with mapQ>=30, then count PCR dups, cis-, trans-,
    args = p.parse_args()
    mh_reads_summary_v2(**vars(args))
