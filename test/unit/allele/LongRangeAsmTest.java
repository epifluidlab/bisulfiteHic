import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import main.java.edu.mit.compbio.bisulfitehic.utils.LongRangeAsm;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

/**
 * End-to-end checks of LongRangeAsm on synthetic data (see Synth.java for the genome layout).
 * SNPs: chr1:10007 and chr2:30007, both A|T (allele A = 'A', allele B = 'T').
 *
 * usage: java LongRangeAsmTest <scratch_dir>      exit status 0 = pass, 1 = fail
 */
public class LongRangeAsmTest {

	static final SAMFileHeader H = Synth.header(SAMFileHeader.SortOrder.coordinate);

	// A1: allele A, methylated, mate 40 kb away.  B1: allele B, unmethylated, 7 leading soft-clipped bases.
	// B2: allele B, unmethylated.  C1: allele A at the chr2 SNP.
	static List<SAMRecord> base(){
		ArrayList<SAMRecord> r = new ArrayList<SAMRecord>();
		Synth.addPair(r, H, "A1", "chr1", 9981, 0, true, 10007, 'A', "chr1", 50001, true, true);
		Synth.addPair(r, H, "B1", "chr1", 9981, 7, false, 10007, 'T', "chr1", 52001, false, true);
		Synth.addPair(r, H, "B2", "chr1", 9981, 0, false, 10007, 'T', "chr1", 54001, false, true);
		Synth.addPair(r, H, "C1", "chr2", 29981, 0, true, 30007, 'A', "chr2", 90001, true, true);
		return r;
	}

	// htsjdk-style TLEN (5'-to-5' distance) instead of the 0 that Bhmem writes for chimeric pairs
	static List<SAMRecord> legacyTlen(List<SAMRecord> recs){
		for(SAMRecord r : recs){
			if(r.getReferenceName().equals(r.getMateReferenceName())){
				r.setInferredInsertSize(r.getMateAlignmentStart() - r.getAlignmentStart());
			}
		}
		return recs;
	}

	static File[] run(File dir, String tag, File vcf, List<SAMRecord> recs, String... extra) throws Exception {
		File bam = new File(dir, tag + ".bam");
		Synth.writeBam(bam, H, recs);
		File tab = new File(dir, tag + ".ASM.tab.gz");
		File ref = new File(dir, tag + ".ref.tab.gz");
		File alt = new File(dir, tag + ".alt.tab.gz");
		ArrayList<String> args = new ArrayList<String>(Arrays.asList(vcf.getPath(), bam.getPath(), tab.getPath(), ref.getPath(), alt.getPath(),
				"-turnOffBisulfiteFilter", "-coverageRef", "1", "-coverageAlt", "1"));
		args.addAll(Arrays.asList(extra));
		new LongRangeAsm().doMain(args.toArray(new String[0]));
		return new File[]{tab, ref, alt};
	}

	static String[] row(File tab, String chr, int pos) throws Exception {
		for(String line : Synth.readGz(tab)){
			if(line.startsWith("#")) continue;
			String[] c = line.split("\t");
			if(c[0].equals(chr) && c[1].equals(String.valueOf(pos))) return c;
		}
		return null;
	}

	// per-CpG file -> "chr:start-end C/total"
	static TreeSet<String> cpgs(File f) throws Exception {
		TreeSet<String> s = new TreeSet<String>();
		for(String line : Synth.readGz(f)){
			String[] c = line.split("\t");
			s.add(c[0] + ":" + c[1] + "-" + c[2] + " " + c[4] + "/" + c[5]);
		}
		return s;
	}

	static TreeSet<String> expect(String chr, int methy, int total, int... cytosines){
		TreeSet<String> s = new TreeSet<String>();
		for(int p : cytosines){
			s.add(chr + ":" + (p - 1) + "-" + p + " " + (methy * total) + "/" + total);
		}
		return s;
	}

	static String get(String[] row, int col1based){
		return row == null || row.length < col1based ? "<missing>" : row[col1based - 1];
	}

	public static void main(String[] argv) throws Exception {
		File dir = Synth.freshDir(argv[0]);
		File vcf = new File(dir, "snps.vcf.gz");
		Synth.writeVcf(vcf, new int[][]{{1, 10007}, {2, 30007}});

		// ---- run 1: BAM with real TLEN, default options ----
		System.out.println("run 1: two chromosomes, soft-clipped read (BAM with real TLEN)");
		try{
			File[] out = run(dir, "run1", vcf, legacyTlen(base()));
			String[] r = row(out[0], "chr1", 10007);
			Synth.check(r != null, "chr1:10007 is reported");
			Synth.check(get(r, 5).equals("1") && get(r, 6).equals("3") && get(r, 7).equals("3"),
					"allele A: 1 pair, 3 methylated of 3 calls (got " + get(r, 5) + "," + get(r, 6) + "," + get(r, 7) + ")");
			Synth.check(get(r, 10).equals("2") && get(r, 11).equals("0") && get(r, 12).equals("5"),
					"allele B: 2 pairs incl. the soft-clipped read, 0 methylated of 5 calls (got " + get(r, 10) + "," + get(r, 11) + "," + get(r, 12) + ")");
			Synth.check(get(r, 15).equals("3") && get(r, 16).equals("3") && get(r, 18).equals("0") && get(r, 19).equals("6"),
					"end 2: allele A 3/3, allele B 0/6");

			List<String> lines = Synth.readGz(out[0]);
			String header = lines.isEmpty() ? "" : lines.get(0);
			Synth.check(r != null && header.split("\t").length == r.length, "header and row have the same number of columns");
			Synth.check(header.contains("AlleleA_CT\t") && !header.contains("AlleleA_T\t"),
					"the column holding C+T totals is not labelled _T");

			TreeSet<String> expRef = expect("chr1", 1, 1, 9983, 10003, 10023, 50003, 50023, 50043);
			expRef.addAll(expect("chr2", 1, 1, 29983, 30003, 30023, 90003, 90023, 90043));
			TreeSet<String> gotRef = cpgs(out[1]);
			Synth.check(gotRef.equals(expRef), "allele A per-CpG file: chr1 CpGs are not repeated under chr2 (got " + gotRef + ")");
			Synth.check(Synth.readGz(out[1]).size() == gotRef.size(), "allele A per-CpG file has no duplicated lines");

			TreeSet<String> expAlt = expect("chr1", 0, 1, 9983, 52003, 52023, 52043, 54003, 54023, 54043);
			expAlt.addAll(expect("chr1", 0, 2, 10003, 10023));
			TreeSet<String> gotAlt = cpgs(out[2]);
			Synth.check(gotAlt.equals(expAlt), "allele B per-CpG file (got " + gotAlt + ")");
		}catch(Exception e){
			Synth.check(false, "run 1 threw " + e);
		}

		// ---- run 2: same pairs but TLEN = 0 as written by Bhmem for chimeric pairs ----
		System.out.println("run 2: long-range cis pairs with TLEN=0, default options");
		try{
			File[] out = run(dir, "run2", vcf, base());
			String[] r = row(out[0], "chr1", 10007);
			Synth.check(r != null && get(r, 5).equals("1") && get(r, 10).equals("2"),
					"pairs 40 kb apart are kept although TLEN is 0 (A=" + get(r, 5) + ", B=" + get(r, 10) + ")");
		}catch(Exception e){
			Synth.check(false, "run 2 threw " + e);
		}

		// ---- run 3: -diffChr (without -useBadMate): short cis pair with TLEN=0 and a trans pair ----
		System.out.println("run 3: -diffChr with a 220 bp cis pair (TLEN=0) and a chr1-chr2 pair");
		try{
			List<SAMRecord> recs = base();
			Synth.addPair(recs, H, "S1", "chr1", 9981, 0, true, 10007, 'A', "chr1", 10201, true, true);
			Synth.addPair(recs, H, "T1", "chr1", 9981, 0, true, 10007, 'A', "chr2", 70001, true, false);
			File[] out = run(dir, "run3", vcf, recs, "-diffChr");
			String[] r = row(out[0], "chr1", 10007);
			Synth.check(get(r, 5).equals("2") && get(r, 6).equals("6") && get(r, 7).equals("6"),
					"allele A = long cis pair + trans pair; the 220 bp pair is excluded (got " + get(r, 5) + "," + get(r, 6) + "," + get(r, 7) + ")");
			TreeSet<String> gotRef = cpgs(out[1]);
			Synth.check(gotRef.containsAll(expect("chr2", 1, 1, 70003, 70023, 70043)),
					"CpGs of the chr2 mate are written under chr2");
			boolean wrongChr = false;
			for(String s : gotRef){
				if(s.startsWith("chr1:700")) wrongChr = true;
			}
			Synth.check(!wrongChr, "CpGs of the chr2 mate are not written under chr1");
			Synth.check(Synth.readGz(out[1]).size() == gotRef.size(), "allele A per-CpG file has no duplicated lines");
		}catch(Exception e){
			Synth.check(false, "run 3 threw " + e);
		}

		// ---- run 4: a read whose mate is absent from the BAM ----
		System.out.println("run 4: mate record missing from the BAM");
		try{
			List<SAMRecord> recs = base();
			recs.add(Synth.read(H, "N1", "chr1", 9981, 50, 0, false, true, 10007, 'A', "chr1", 60001, true));
			File[] out = run(dir, "run4", vcf, legacyTlen(recs));
			String[] r = row(out[0], "chr1", 10007);
			Synth.check(get(r, 5).equals("1") && get(r, 22).equals("2") && get(r, 23).equals("6") && get(r, 24).equals("6"),
					"read without mate counts only in the local columns (A=" + get(r, 5) + ", A_local=" + get(r, 22) + ")");
		}catch(Exception e){
			Synth.check(false, "run 4 threw " + e);
		}

		Synth.finish("LongRangeAsmTest");
	}
}
