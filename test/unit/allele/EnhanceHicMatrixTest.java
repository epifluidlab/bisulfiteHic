import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import main.java.edu.mit.compbio.bisulfitehic.utils.EnhanceHicMatrix;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * End-to-end checks of EnhanceHicMatrix -bedpe on synthetic data (see Synth.java).
 * Four read pairs link chr1:12001-12350 with chr1:62001-62350 and overlap the CpG list;
 * a fifth pair (14001 / 64001) lies in the same anchors but overlaps no listed CpG.
 *
 * usage: java EnhanceHicMatrixTest <scratch_dir>      exit status 0 = pass, 1 = fail
 */
public class EnhanceHicMatrixTest {

	static List<String> run(File dir, String tag, File sizes, File cpg, File bam, String bedpeLine) throws Exception {
		File bedpe = new File(dir, tag + ".bedpe");
		Synth.writeText(bedpe, bedpeLine + "\n");
		String prefix = new File(dir, tag).getPath();
		new EnhanceHicMatrix().doMain(new String[]{prefix, sizes.getPath(), cpg.getPath(), bam.getPath(), "-bedpe", bedpe.getPath(), "-turnOffBisulfiteFilter"});
		return Synth.readGz(new File(prefix + ".enhancedHicMatrix.txt.gz"));
	}

	public static void main(String[] argv) throws Exception {
		File dir = Synth.freshDir(argv[0]);
		SAMFileHeader h = Synth.header(SAMFileHeader.SortOrder.coordinate);
		ArrayList<SAMRecord> recs = new ArrayList<SAMRecord>();
		boolean[] m1 = {true, false, true, false};
		boolean[] m2 = {true, false, true, true};
		for(int i = 0; i < 4; i++){
			Synth.addPair(recs, h, "p" + i, "chr1", 12001 + 100 * i, 0, m1[i], 0, 'N', "chr1", 62001 + 100 * i, m2[i], true);
		}
		Synth.addPair(recs, h, "outside", "chr1", 14001, 0, false, 0, 'N', "chr1", 64001, false, true);
		File bam = new File(dir, "reads.bam");
		Synth.writeBam(bam, h, recs);
		File sizes = new File(dir, "chrom.sizes");
		Synth.writeText(sizes, "chr1\t" + Synth.CONTIG_LEN + "\nchr2\t" + Synth.CONTIG_LEN + "\n");
		File cpg = new File(dir, "cpg_list.bed");
		Synth.writeText(cpg, "chr1\t12000\t12400\nchr1\t62000\t62400\n");

		System.out.println("run 1: bedpe anchors exactly as wide as -window (5 kb)");
		try{
			List<String> lines = run(dir, "narrow", sizes, cpg, bam, "chr1\t10000\t15000\tchr1\t60000\t65000");
			Synth.check(lines.size() == 1, "an anchor pair of width == window is analysed (output lines: " + lines.size() + ")");
		}catch(Exception e){
			Synth.check(false, "run 1 threw " + e);
		}

		System.out.println("run 2: 10 kb anchors, one read pair overlaps no CpG of cpg_list.bed");
		try{
			List<String> lines = run(dir, "wide", sizes, cpg, bam, "chr1\t10000\t20000\tchr1\t60000\t70000");
			String n = lines.isEmpty() ? "<no output>" : lines.get(0).split("\t")[9];
			Synth.check(lines.size() == 1 && Double.parseDouble(n) == 4.0, "reads without a listed CpG are not used: 4 read pairs (got " + n + ")");
		}catch(Exception e){
			Synth.check(false, "run 2 threw " + e);
		}

		Synth.finish("EnhanceHicMatrixTest");
	}
}
