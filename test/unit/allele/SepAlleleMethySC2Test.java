import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import main.java.edu.mit.compbio.bisulfitehic.utils.SepAlleleMethySC2;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

/**
 * End-to-end checks of SepAlleleMethySC2 on a synthetic read-name-sorted BAM (see Synth.java).
 * SNP chr1:10007 A|T.  Pair p1 = allele A, methylated; pair p2 = allele B, unmethylated.
 * The SNP sits at read offset 26 of end 1; end 2 (at 10097) has a CpG at its own offset 26.
 *
 * usage: java SepAlleleMethySC2Test <scratch_dir>      exit status 0 = pass, 1 = fail
 */
public class SepAlleleMethySC2Test {

	static final SAMFileHeader H = Synth.header(SAMFileHeader.SortOrder.queryname);

	static List<SAMRecord> pairs(){
		ArrayList<SAMRecord> r = new ArrayList<SAMRecord>();
		Synth.addPair(r, H, "p1", "chr1", 9981, 0, true, 10007, 'A', "chr1", 10097, true, true);
		Synth.addPair(r, H, "p2", "chr1", 9981, 0, false, 10007, 'T', "chr1", 10097, false, true);
		return r;
	}

	static TreeSet<String> expect(int methy){
		TreeSet<String> s = new TreeSet<String>();
		for(int p : new int[]{9983, 10003, 10023, 10103, 10123, 10143}){
			s.add("chr1\t" + (p - 1) + "\t" + p + "\t" + methy);
		}
		return s;
	}

	static void runAndCheck(File dir, String tag, File vcf, List<SAMRecord> recs, String what){
		try{
			File bam = new File(dir, tag + ".bam");
			Synth.writeBam(bam, H, recs);
			String prefix = new File(dir, tag).getPath();
			new SepAlleleMethySC2().doMain(new String[]{prefix, vcf.getPath(), bam.getPath(), "-turnOffBisulfiteFilter"});
			TreeSet<String> a = new TreeSet<String>(Synth.readGz(new File(prefix + ".alleleA_methy.bedgraph.gz")));
			TreeSet<String> b = new TreeSet<String>(Synth.readGz(new File(prefix + ".alleleB_methy.bedgraph.gz")));
			Synth.check(a.equals(expect(1)), what + ": allele A CpGs (got " + a.toString().replace('\t', ' ') + ")");
			Synth.check(b.equals(expect(0)), what + ": allele B CpGs (got " + b.toString().replace('\t', ' ') + ")");
		}catch(Exception e){
			Synth.check(false, what + ": threw " + e);
		}
	}

	public static void main(String[] argv) throws Exception {
		File dir = Synth.freshDir(argv[0]);
		File vcf = new File(dir, "snps.vcf.gz");
		Synth.writeVcf(vcf, new int[][]{{1, 10007}});

		System.out.println("run 1: two clean pairs");
		runAndCheck(dir, "clean", vcf, pairs(), "mate CpG at the SNP's read offset is kept");

		System.out.println("run 2: a read without its mate precedes the pairs");
		ArrayList<SAMRecord> recs = new ArrayList<SAMRecord>();
		recs.add(Synth.read(H, "a0", "chr1", 20001, 50, 0, false, true, 0, 'N', "chr1", 60001, true));
		recs.addAll(pairs());
		runAndCheck(dir, "singleton", vcf, recs, "mates are matched by read name");

		// end 1 is 16S34M starting at 9997: SNP at offset 26 of the whole read, CpG 10023 at offset 26 from the first aligned base
		System.out.println("run 3: soft-clipped end 1");
		try{
			ArrayList<SAMRecord> clip = new ArrayList<SAMRecord>();
			Synth.addPair(clip, H, "p3", "chr1", 9981, 16, true, 10007, 'A', "chr1", 10097, true, true);
			File bam = new File(dir, "clip.bam");
			Synth.writeBam(bam, H, clip);
			String prefix = new File(dir, "clip").getPath();
			new SepAlleleMethySC2().doMain(new String[]{prefix, vcf.getPath(), bam.getPath(), "-turnOffBisulfiteFilter"});
			TreeSet<String> a = new TreeSet<String>(Synth.readGz(new File(prefix + ".alleleA_methy.bedgraph.gz")));
			Synth.check(a.contains("chr1\t10022\t10023\t1") && a.contains("chr1\t10002\t10003\t1") && a.size() == 5,
					"all CpGs of a soft-clipped read are reported (got " + a.toString().replace('\t', ' ') + ")");
		}catch(Exception e){
			Synth.check(false, "run 3 threw " + e);
		}

		Synth.finish("SepAlleleMethySC2Test");
	}
}
