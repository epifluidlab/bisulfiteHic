import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Tiny synthetic-data builder shared by the tests in this folder.
 *
 * The "genome" is two contigs (chr1, chr2; 200 kb each) whose sequence is PATTERN repeated,
 * so every reference base is known without a FASTA file:
 *   - one CpG per 20 bp, its C at 1-based positions p with (p-1)%20 == 2   (3, 23, 43, ...)
 *   - an 'A' that is used as an A/T SNP site at positions p with (p-1)%20 == 6 (7, 27, ...)
 * All reads are built on the C->T bisulfite strand (first-of-pair forward, or second-of-pair
 * reverse), carry an MD tag, and have TLEN = 0, which is what Bhmem writes for chimeric
 * Hi-C pairs.
 */
public class Synth {

	public static final String PATTERN = "AACGTTAGTCATGAATTGAA";
	public static final int CONTIG_LEN = 200000;

	public static int failures = 0;

	public static void check(boolean ok, String what){
		if(ok){
			System.out.println("  ok    " + what);
		}else{
			failures++;
			System.out.println("  FAIL  " + what);
		}
	}

	public static void finish(String testName){
		if(failures > 0){
			System.out.println(testName + ": " + failures + " check(s) FAILED");
			System.exit(1);
		}
		System.out.println(testName + ": all checks passed");
		System.exit(0);
	}

	public static char ref(int pos){
		return PATTERN.charAt((pos - 1) % PATTERN.length());
	}

	public static SAMSequenceDictionary dict(){
		SAMSequenceDictionary d = new SAMSequenceDictionary();
		d.addSequence(new SAMSequenceRecord("chr1", CONTIG_LEN));
		d.addSequence(new SAMSequenceRecord("chr2", CONTIG_LEN));
		return d;
	}

	public static SAMFileHeader header(SAMFileHeader.SortOrder order){
		SAMFileHeader h = new SAMFileHeader();
		h.setSequenceDictionary(dict());
		h.setSortOrder(order);
		return h;
	}

	/**
	 * A read on the C->T strand.
	 * @param softClip   number of leading soft-clipped bases (CIGAR softClip S len M)
	 * @param second     second of pair; the strand flag is set so the read is on the C->T strand
	 * @param methylated CpG cytosines read as C (true) or T (false); non-CpG C always reads as T
	 * @param snpPos     1-based reference position whose base is replaced by snpBase (0 = none)
	 * @param proper     proper-pair flag
	 */
	public static SAMRecord read(SAMFileHeader h, String name, String chr, int start, int len, int softClip, boolean second,
			boolean methylated, int snpPos, char snpBase, String mateChr, int mateStart, boolean proper){
		StringBuilder seq = new StringBuilder();
		for(int i = 0; i < softClip; i++){
			seq.append('G');
		}
		StringBuilder md = new StringBuilder();
		int match = 0;
		for(int pos = start; pos < start + len; pos++){
			char r = ref(pos);
			char b = r;
			if(pos == snpPos){
				b = snpBase;
			}else if(r == 'C'){
				b = (ref(pos + 1) == 'G' && methylated) ? 'C' : 'T';
			}
			seq.append(b);
			if(b == r){
				match++;
			}else{
				md.append(match).append(r);
				match = 0;
			}
		}
		md.append(match);
		char[] q = new char[seq.length()];
		Arrays.fill(q, 'I');

		SAMRecord rec = new SAMRecord(h);
		rec.setReadName(name);
		rec.setReferenceName(chr);
		rec.setAlignmentStart(start);
		rec.setCigarString(softClip > 0 ? (softClip + "S" + len + "M") : (len + "M"));
		rec.setReadString(seq.toString());
		rec.setBaseQualityString(new String(q));
		rec.setMappingQuality(60);
		rec.setReadPairedFlag(true);
		rec.setFirstOfPairFlag(!second);
		rec.setSecondOfPairFlag(second);
		rec.setReadNegativeStrandFlag(second);
		rec.setMateNegativeStrandFlag(!second);
		rec.setProperPairFlag(proper);
		rec.setMateUnmappedFlag(false);
		rec.setMateReferenceName(mateChr);
		rec.setMateAlignmentStart(mateStart);
		rec.setInferredInsertSize(0);
		rec.setAttribute("MD", md.toString());
		return rec;
	}

	/** Both ends of a pair: end 1 carries the SNP allele (if snpPos > 0), end 2 is the distal mate. */
	public static void addPair(List<SAMRecord> out, SAMFileHeader h, String name, String chr1, int start1, int softClip1, boolean meth1,
			int snpPos, char snpBase, String chr2, int start2, boolean meth2, boolean proper){
		int len1 = 50 - softClip1;
		int alnStart1 = start1 + softClip1;
		out.add(read(h, name, chr1, alnStart1, len1, softClip1, false, meth1, snpPos, snpBase, chr2, start2, proper));
		out.add(read(h, name, chr2, start2, 50, 0, true, meth2, 0, 'N', chr1, alnStart1, proper));
	}

	public static void writeBam(File f, SAMFileHeader h, List<SAMRecord> recs){
		boolean coordinate = h.getSortOrder() == SAMFileHeader.SortOrder.coordinate;
		SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(coordinate).makeBAMWriter(h, !coordinate, f);
		for(SAMRecord r : recs){
			w.addAlignment(r);
		}
		w.close();
	}

	/** Phased het A|T SNPs (genotype 0|1: allele A = reference 'A', allele B = 'T'); each entry is {contigIndex(1 or 2), pos}. */
	public static void writeVcf(File f, int[][] snps){
		SAMSequenceDictionary d = dict();
		HashSet<VCFHeaderLine> meta = new HashSet<VCFHeaderLine>();
		meta.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
		VCFHeader vh = new VCFHeader(meta, Arrays.asList("S1"));
		vh.setSequenceDictionary(d);
		VariantContextWriter w = new VariantContextWriterBuilder().setOutputFile(f).setReferenceDictionary(d).build();
		w.writeHeader(vh);
		for(int[] s : snps){
			Allele a = Allele.create("A", true);
			Allele t = Allele.create("T", false);
			w.add(new VariantContextBuilder("synth", "chr" + s[0], s[1], s[1], Arrays.asList(a, t))
					.genotypes(new GenotypeBuilder("S1", Arrays.asList(a, t)).phased(true).make()).make());
		}
		w.close();
	}

	public static void writeText(File f, String text) throws Exception {
		FileWriter w = new FileWriter(f);
		w.write(text);
		w.close();
	}

	public static List<String> readGz(File f) throws Exception {
		ArrayList<String> lines = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f))));
		String line;
		while((line = br.readLine()) != null){
			lines.add(line);
		}
		br.close();
		return lines;
	}

	public static File freshDir(String path){
		File d = new File(path);
		d.mkdirs();
		File[] old = d.listFiles();
		if(old != null){
			for(File o : old){
				if(o.isFile()){
					o.delete();
				}
			}
		}
		return d;
	}
}
