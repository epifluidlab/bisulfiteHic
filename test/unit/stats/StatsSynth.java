import java.io.File;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

/**
 * Synthetic bisulfite read pairs for the BAM-based tests in this folder.
 *
 * The reference is 'A' everywhere except a CG at every listed CpG position (C at the 1-based
 * position, G right after it). Mate 1 is first-of-pair on the forward strand and mate 2 is
 * second-of-pair on the reverse strand, so both are C->T reads: a methylated C reads C, an
 * unmethylated C reads T. The reference is only carried by the MD tag (the tools rebuild it
 * from there), so no FASTA is needed.
 */
public class StatsSynth {
	static final String CHR = "chr1";
	static final int CHR_LEN = 100000;

	static SAMFileHeader header(SAMFileHeader.SortOrder order) {
		SAMFileHeader h = new SAMFileHeader();
		h.addSequence(new SAMSequenceRecord(CHR, CHR_LEN));
		h.setSortOrder(order);
		return h;
	}

	static char refBase(int pos, int[] cpgs) {
		for (int c : cpgs) {
			if (pos == c) {
				return 'C';
			}
			if (pos == c + 1) {
				return 'G';
			}
		}
		return 'A';
	}

	/** base a C->T read shows at pos: T for an unmethylated CpG cytosine, snpBase at snpPos, else the reference */
	static char readBase(int pos, int[] cpgs, boolean[] meth, int snpPos, char snpBase) {
		if (pos == snpPos) {
			return snpBase;
		}
		for (int k = 0; k < cpgs.length; k++) {
			if (pos == cpgs[k] && !meth[k]) {
				return 'T';
			}
		}
		return refBase(pos, cpgs);
	}

	/**
	 * One mate. The alignment is lenBefore matched bases, then (if delLen > 0) a deletion of delLen
	 * reference bases, then lenAfter matched bases. meth[k] is the state of cpgs[k]; snpPos 0 = no SNP.
	 */
	static SAMRecord mate(SAMFileHeader h, String name, boolean first, int start, int lenBefore, int delLen, int lenAfter,
			int[] cpgs, boolean[] meth, int snpPos, char snpBase, int mateStart) {
		StringBuilder seq = new StringBuilder();
		StringBuilder md = new StringBuilder();
		int run = 0;
		int pos = start;
		for (int block = 0; block < 2; block++) {
			int len = block == 0 ? lenBefore : lenAfter;
			for (int i = 0; i < len; i++, pos++) {
				char ref = refBase(pos, cpgs);
				char base = readBase(pos, cpgs, meth, snpPos, snpBase);
				seq.append(base);
				if (base == ref) {
					run++;
				} else {
					md.append(run).append(ref);
					run = 0;
				}
			}
			if (block == 0 && delLen > 0) {
				md.append(run).append('^');
				for (int i = 0; i < delLen; i++, pos++) {
					md.append(refBase(pos, cpgs));
				}
				run = 0;
			}
		}
		md.append(run);
		String cigar = delLen > 0 ? lenBefore + "M" + delLen + "D" + lenAfter + "M" : (lenBefore + lenAfter) + "M";
		char[] qual = new char[seq.length()];
		Arrays.fill(qual, 'I');

		SAMRecord r = new SAMRecord(h);
		r.setReadName(name);
		r.setReferenceName(CHR);
		r.setAlignmentStart(start);
		r.setCigarString(cigar);
		r.setReadString(seq.toString());
		r.setBaseQualityString(new String(qual));
		r.setMappingQuality(60);
		r.setReadPairedFlag(true);
		r.setProperPairFlag(true);
		r.setFirstOfPairFlag(first);
		r.setSecondOfPairFlag(!first);
		r.setReadNegativeStrandFlag(!first);
		r.setMateNegativeStrandFlag(first);
		r.setMateReferenceName(CHR);
		r.setMateAlignmentStart(mateStart);
		r.setAttribute("MD", md.toString());
		return r;
	}

	/** ungapped 40 bp mate that starts 10 bp upstream of anchor */
	static SAMRecord mate(SAMFileHeader h, String name, boolean first, int anchor, int[] cpgs, boolean[] meth, int snpPos, char snpBase, int mateAnchor) {
		return mate(h, name, first, anchor - 10, 40, 0, 0, cpgs, meth, snpPos, snpBase, mateAnchor - 10);
	}

	/** coordinate-sorted BAM with a .bai index */
	static void writeBam(File f, SAMFileHeader h, List<SAMRecord> recs) {
		SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(h, false, f);
		for (SAMRecord r : recs) {
			w.addAlignment(r);
		}
		w.close();
	}

	/** SAM text file that keeps the given record order (the two mates of a pair next to each other) */
	static void writeSam(File f, SAMFileHeader h, List<SAMRecord> recs) {
		SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(h, true, f);
		for (SAMRecord r : recs) {
			w.addAlignment(r);
		}
		w.close();
	}

	/** bgzipped, tabix-indexed VCF (f must end with .vcf.gz) with one heterozygous SNP ref/alt at pos for sample S1 */
	static void writeHetSnpVcf(File f, int pos, String ref, String alt) {
		SAMFileHeader h = header(SAMFileHeader.SortOrder.coordinate);
		VCFHeader vh = new VCFHeader(new LinkedHashSet<VCFHeaderLine>(), Arrays.asList("S1"));
		vh.addMetaDataLine(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
		vh.setSequenceDictionary(h.getSequenceDictionary());
		VariantContextWriter w = new VariantContextWriterBuilder().setOutputFile(f)
				.setReferenceDictionary(h.getSequenceDictionary()).setOption(Options.INDEX_ON_THE_FLY).build();
		w.writeHeader(vh);
		List<Allele> alleles = Arrays.asList(Allele.create(ref, true), Allele.create(alt, false));
		w.add(new VariantContextBuilder("synth", CHR, pos, pos, alleles).genotypes(new GenotypeBuilder("S1", alleles).make()).make());
		w.close();
	}
}
