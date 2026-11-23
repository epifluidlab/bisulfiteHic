import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import main.java.edu.mit.compbio.bisulfitehic.utils.SepAlleleMethySC;

import java.io.File;
import java.util.ArrayList;
import java.util.TreeSet;

/**
 * End-to-end check of SepAlleleMethySC on synthetic data (see Synth.java).
 * SNP chr1:10007 A|T.  One allele-A read, methylated, CIGAR 16S34M starting at 9997:
 * the SNP is at offset 26 of the whole read (10 from the first aligned base) and the
 * CpG at 10023 is at offset 26 from the first aligned base.
 *
 * usage: java SepAlleleMethySCTest <scratch_dir>      exit status 0 = pass, 1 = fail
 */
public class SepAlleleMethySCTest {

	public static void main(String[] argv) throws Exception {
		File dir = Synth.freshDir(argv[0]);
		SAMFileHeader h = Synth.header(SAMFileHeader.SortOrder.coordinate);
		File vcf = new File(dir, "snps.vcf.gz");
		Synth.writeVcf(vcf, new int[][]{{1, 10007}});
		ArrayList<SAMRecord> recs = new ArrayList<SAMRecord>();
		Synth.addPair(recs, h, "c1", "chr1", 9981, 16, true, 10007, 'A', "chr1", 10097, true, true);
		File bam = new File(dir, "reads.bam");
		Synth.writeBam(bam, h, recs);

		System.out.println("run 1: soft-clipped read, a CpG sits at (SNP read offset) bases after the first aligned base");
		try{
			String prefix = new File(dir, "out").getPath();
			new SepAlleleMethySC().doMain(new String[]{prefix, vcf.getPath(), bam.getPath(), "-turnOffBisulfiteFilter"});
			TreeSet<String> a = new TreeSet<String>(Synth.readGz(new File(prefix + ".alleleA_methy.bedgraph.gz")));
			TreeSet<String> exp = new TreeSet<String>();
			exp.add("chr1\t10002\t10003\t1");
			exp.add("chr1\t10022\t10023\t1");
			Synth.check(a.equals(exp), "both CpGs of the soft-clipped read are reported (got " + a.toString().replace('\t', ' ') + ")");
		}catch(Exception e){
			Synth.check(false, "run 1 threw " + e);
		}

		Synth.finish("SepAlleleMethySCTest");
	}
}
