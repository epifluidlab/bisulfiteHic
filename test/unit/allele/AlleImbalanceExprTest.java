import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import main.java.edu.mit.compbio.bisulfitehic.utils.AlleImbalanceExpr;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * End-to-end checks of AlleImbalanceExpr on synthetic data (see Synth.java for the genome layout).
 * SNPs: chr1:10007 (6 reads allele A, 1 soft-clipped read allele B) and chr1:12007 (30 reads allele A).
 *
 * usage: java AlleImbalanceExprTest <scratch_dir>      exit status 0 = pass, 1 = fail
 */
public class AlleImbalanceExprTest {

	static List<String[]> rows(File f, String gene) throws Exception {
		ArrayList<String[]> r = new ArrayList<String[]>();
		for(String line : Synth.readGz(f)){
			if(line.startsWith("#")) continue;
			String[] c = line.split("\t");
			if(c[3].equals(gene)) r.add(c);
		}
		return r;
	}

	public static void main(String[] argv) throws Exception {
		File dir = Synth.freshDir(argv[0]);
		SAMFileHeader h = Synth.header(SAMFileHeader.SortOrder.coordinate);
		File vcf = new File(dir, "snps.vcf.gz");
		Synth.writeVcf(vcf, new int[][]{{1, 10007}, {1, 12007}});

		ArrayList<SAMRecord> recs = new ArrayList<SAMRecord>();
		for(int i = 0; i < 6; i++){
			recs.add(Synth.read(h, "a" + i, "chr1", 9981, 50, 0, false, true, 10007, 'A', "chr1", 50001, true));
		}
		recs.add(Synth.read(h, "b0", "chr1", 9988, 43, 7, false, true, 10007, 'T', "chr1", 50001, true));
		for(int i = 0; i < 30; i++){
			recs.add(Synth.read(h, "z" + i, "chr1", 11981, 50, 0, false, true, 12007, 'A', "chr1", 50001, true));
		}
		File bam = new File(dir, "reads.bam");
		Synth.writeBam(bam, h, recs);

		// geneX and geneY share coordinates; geneEdge starts (BED, 0-based) right after the SNP at 10007
		File bed = new File(dir, "genes.bed");
		Synth.writeText(bed, "chr1\t9000\t11000\tgeneX\t.\t+\n"
				+ "chr1\t9000\t11000\tgeneY\t.\t-\n"
				+ "chr1\t11500\t12500\tgeneZ\t.\t+\n"
				+ "chr1\t10007\t10100\tgeneEdge\t.\t+\n");

		System.out.println("run 1: one line per SNP");
		try{
			File out = new File(dir, "perSnp.tab.gz");
			new AlleImbalanceExpr().doMain(new String[]{bed.getPath(), vcf.getPath(), bam.getPath(), out.getPath()});
			List<String[]> x = rows(out, "geneX");
			Synth.check(x.size() == 1, "geneX is reported once");
			String a = x.isEmpty() ? "<missing>" : x.get(0)[10];
			String b = x.isEmpty() ? "<missing>" : x.get(0)[12];
			Synth.check(a.equals("6") && b.equals("1"), "soft-clipped read is counted for its own allele: A=6, B=1 (got A=" + a + ", B=" + b + ")");
			Synth.check(rows(out, "geneY").size() == 1, "geneY (same coordinates as geneX) is reported too");
			Synth.check(rows(out, "geneEdge").isEmpty(), "SNP one base before a BED interval is not assigned to it");
			List<String[]> z = rows(out, "geneZ");
			double p = z.isEmpty() ? Double.NaN : Double.parseDouble(z.get(0)[13]);
			Synth.check(p > 0 && p < 1e-6, "p-value of a 30:0 SNP is not rounded to 0 (got " + (z.isEmpty() ? "<missing>" : z.get(0)[13]) + ")");
			String header = Synth.readGz(out).get(0);
			Synth.check(!x.isEmpty() && header.split("\t").length == x.get(0).length, "header and rows have the same number of columns");
		}catch(Exception e){
			Synth.check(false, "run 1 threw " + e);
		}

		System.out.println("run 2: -mergedInSameGene");
		try{
			File out = new File(dir, "perGene.tab.gz");
			new AlleImbalanceExpr().doMain(new String[]{bed.getPath(), vcf.getPath(), bam.getPath(), out.getPath(), "-mergedInSameGene"});
			List<String[]> x = rows(out, "geneX");
			String[] header = Synth.readGz(out).get(0).split("\t");
			Synth.check(!x.isEmpty() && header.length == x.get(0).length,
					"header and rows have the same number of columns (header " + header.length + ", row " + (x.isEmpty() ? 0 : x.get(0).length) + ")");
			int ia = -1;
			int ib = -1;
			for(int i = 0; i < header.length; i++){
				if(header[i].equals("AlleleA_count")) ia = i;
				if(header[i].equals("AlleleB_count")) ib = i;
			}
			Synth.check(!x.isEmpty() && ia >= 0 && ib >= 0 && ib < x.get(0).length && x.get(0)[ia].equals("6") && x.get(0)[ib].equals("1"),
					"columns named AlleleA_count/AlleleB_count hold 6 and 1");
		}catch(Exception e){
			Synth.check(false, "run 2 threw " + e);
		}

		Synth.finish("AlleImbalanceExprTest");
	}
}
