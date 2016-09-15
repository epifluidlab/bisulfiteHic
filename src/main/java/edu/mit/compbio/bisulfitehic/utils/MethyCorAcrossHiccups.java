/**
 * MethyCorAcrossHiccups.java
 * Sep 12, 2016
 * 1:03:54 PM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import main.java.edu.mit.compbio.bissnp2.utils.BisSnp2Utils;
import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.biojava.nbio.genome.parsers.twobit.TwoBitParser;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * accept BEDPE file. for each rows, locate reads from paired end bam file within the region, also locate SNP from VCF file within the region
 * for region within 20kb, skip it
 * iterate each reads that two ends are located at two regions, each fragment is a datapoint in the methy_cor output txt.gz file. (chr1,read.start1,read.end1,chr2,read.start2,read.end2,region_name, id, hiccups_ob, hiccups_exp, hiccups_fdr, methy_end1, cpg_num_end1, methy_end2, cpg_num_end2) 
 * iterate each reads that overlapped with heaterozygous SNPs in VCF file (also have another end mapped to the other part of hiccups.) summaryize the methy ratio on two different allele in two end
 * SNP.chr, SNP.start,SNP.end,region_name, id,hiccups_ob, hiccups_exp, hiccups_fdr, SNP_end.alleleA_all.methy_count,SNP_end.alleleA_all.unmethy_count,SNP_end.alleleB_all.methy_count,SNP_end.alleleB_all.unmethy_count,methyDiff(A-B),fisher_test_p_value, FDR, 
 * SNP_end.alleleA.methy_count,SNP_end.alleleA.unmethy_count,SNP_end.alleleB.methy_count,SNP_end.alleleB.unmethy_count,methyDiff(A-B),fisher_test_p_value, FDR,
 * nonSNP_end.alleleA.methy_count,nonSNP_end.alleleA.unmethy_count,nonSNP_end.alleleB.methy_count,nonSNP_end.alleleB.unmethy_count,methyDiff(A-B), fisher_test_p_value, FDR
 * (maybe further consider the other end have SNP also or not...)
 * then overlapped with GWAS other SNP loci...
 * currently only work for hiccups within the same chromosomes
 */

public class MethyCorAcrossHiccups {

	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;
	
	@Option(name="-turnOffBisulfiteFilter",usage="turn off bisulfite conversion reads filter. Default: false")
	public boolean turnOffBisulfiteFilter = false;
	
	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = -1;

	@Option(name="-minDistEnds",usage="minimum distance between two end to be considered as HiC links. Default: 20000")
	public int minDistEnds = 20000;
	
	@Option(name="-vcf",usage="vcf file. Default: null")
	public String vcf = null;
	
	@Option(name="-coverageRef",usage="For ASM mode, specify the minimum  coverage in SNP's allele Ref, default: 2")
	public int coverageRef = 2;

	@Option(name="-coverageAlt",usage="For ASM mode, specify the minimum coverage in SNP's allele Alt, default: 2")
	public int coverageAlt = 2;
	
	@Option(name="-sample",usage="specify the sample number in the VCF file if contain multiple sample genotype. number started at 1. default: 1")
	public int sample = 1;

	@Option(name="-randomHiccupsInput",usage="when the input file is random hiccups links for the control, no need to restrict each end of fragment at the same hiccups links.., default: not enabled")
	public boolean randomHiccupsInput = false;

	@Option(name="-onlyShort",usage="only use the paired end within very short region. when enable it with positive number, it will automate disable -minDistEnds option, default: -1")
	public int onlyShort = -1;

	@Option(name="-randomSampleReads",usage="only use part of the reads in each interval. good for -randomHiccupsInput mode or -onlyShort mode, default: -1. not enabled")
	public int randomSampleReads = -1;
	
	@Option(name="-skipProcessFirstRow",usage="don't process first row, default: not enabled")
	public boolean skipProcessFirstRow = false;
	
	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "MethyCorAcrossHiccups [opts] outputPrefix mm9.2bit hiccups.bedpe.bed[.gz] input.bam";
	
	private static Logger log = Logger.getLogger(MethyCorAcrossHiccups.class);

	private OutputStreamWriter methyCorWriter = null; 
	private OutputStreamWriter asmWriter = null; 
	private TwoBitParser refParser = null;
	private String preChr="";
	private MersenneTwister generator = null;
	
	private VCFFileReader vcfReader = null;
	
	private static long startTime = -1;
	private static long lineNum=0;
	private static long queryNum=0;
	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
		MethyCorAcrossHiccups mcah = new MethyCorAcrossHiccups();
			BasicConfigurator.configure();
			mcah.doMain(args);
		}
		

		public void doMain(String[] args)
				throws Exception {

						CmdLineParser parser = new CmdLineParser(this);
						//parser.setUsageWidth(80);
						try
						{
							if(help || args.length < 4) throw new CmdLineException(parser, USAGE, new Throwable());
							parser.parseArgument(args);
							
						
						}
						catch (CmdLineException e)
						{
							System.err.println(e.getMessage());
							// print the list of available options
							parser.printUsage(System.err);
							System.err.println();
							return;
						}

						//read input bed file, for each row,
						//String intervalFile = arguments.get(0);
						String outputPrefix = arguments.get(0);
						String refFile = arguments.get(1);
						String hiccupsFile = arguments.get(2);
						String bamFile = arguments.get(3);

						initiate(outputPrefix, refFile);
						SamReader wgsReader = SamReaderFactory.makeDefault().setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true).validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
						log.info("Parsing input hiccups links BEDPE format bed file ...");
						GZIPInputStream gzipInputStream = null;
						BufferedReader br;
						if(hiccupsFile.endsWith(".gz")){
							gzipInputStream = new GZIPInputStream(new FileInputStream(hiccupsFile));
							br = new BufferedReader(new InputStreamReader(gzipInputStream));
							
						}else{
							br = new BufferedReader(new FileReader(hiccupsFile));
						}
						String line;
						
						while( (line = br.readLine()) != null){
							if((skipProcessFirstRow && lineNum==0) || line.startsWith("#")){
								continue;
							}else{
								String[] splitin = line.split("\t");
								String chr1 = splitin[0];
								int start1 = Integer.parseInt(splitin[1]);
								int end1 = Integer.parseInt(splitin[2]);
								String chr2 = splitin[3];
								int start2 = Integer.parseInt(splitin[4]);
								int end2 = Integer.parseInt(splitin[5]);
								if((!chr1.equalsIgnoreCase(chr2)) || (onlyShort<0 && Math.abs(start1-start2)<minDistEnds) || (onlyShort<0 && Math.abs(end1-end2)<minDistEnds)){
									continue;
								}
								double ob = Double.parseDouble(splitin[7]);
								double exp = Double.parseDouble(splitin[11]);
								double fdr = Double.parseDouble(splitin[15]);
								
								if(vcf != null){
									asmCall(chr1, start1, end1, start2, end2, wgsReader, ob, exp, fdr);
								}else{
									if(randomHiccupsInput){
										//log.info("randomHiccupsInput mode: ");
										methyCorrelationRandom(chr1, start1, end1, start2, end2, wgsReader, ob, exp, fdr, bamFile);
									}else{
										if(onlyShort>0){
											//log.info("only use paired end reads within short distance: ");
											methyCorrelationShort(chr1, start1, end1, start2, end2, wgsReader, ob, exp, fdr, bamFile);
										}else{
											methyCorrelation(chr1, start1, end1, start2, end2, wgsReader, ob, exp, fdr, bamFile);
										}
										
									}
								}
								
								
								
							}
								
							
							lineNum++;
							if(lineNum % 1000 == 0){
								log.info("Processing line: " + lineNum);
								
							}
							
						}
						if(hiccupsFile.endsWith(".gz")){
							gzipInputStream.close();
						}
						br.close();
						
						finish();
		}
		
		//currently, not consider condition that each end overlapped a SNP
		//TODO
		private void asmCall(String chr, int start1, int end1, int start2, int end2, SamReader wgsReader, double ob, double exp, double fdr){ 
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + start2 + ":" + end2;
			
			//read SNP list
			//region1
			CloseableIterator<VariantContext> itSnp = vcfReader.query(chr, start1, end1);
			SAMSequenceDictionary seqDict = VCFFileReader.getSequenceDictionary(new File(vcf));
			GenomeLocParser glp = new GenomeLocParser(seqDict);
			HashMap<GenomeLoc, Pair<Byte, Byte>> genotypeMapList1 = new HashMap<GenomeLoc, Pair<Byte, Byte>>();
			while(itSnp.hasNext()){
				VariantContext vc = itSnp.next();
				if(vc.isFiltered() || vc.isIndel() || vc.isMixed() || !vc.isPointEvent())
					continue;
				
					Genotype gt = vc.getGenotype(sample);
					
					if(!gt.isFiltered() && gt.isHet()){ //only phased Het SNP is considered
						GenomeLoc genomeLoc = glp.createGenomeLoc(vc.getContig(),vc.getStart(), vc.getEnd());
						Byte alleleA = gt.getAllele(0).getBases()[0];
						Byte alleleB = gt.getAllele(1).getBases()[0];
						Pair<Byte, Byte> tmp = new Pair<Byte, Byte>(alleleA, alleleB);
						genotypeMapList1.put(genomeLoc, tmp);
					}
				
				
			}
			itSnp.close();
			
			//region2
			itSnp = vcfReader.query(chr, start2, end2);
			HashMap<GenomeLoc, Pair<Byte, Byte>> genotypeMapList2 = new HashMap<GenomeLoc, Pair<Byte, Byte>>();
			while(itSnp.hasNext()){
				VariantContext vc = itSnp.next();
				if(vc.isFiltered() || vc.isIndel() || vc.isMixed() || !vc.isPointEvent())
					continue;
					Genotype gt = vc.getGenotype(sample);				
					if(!gt.isFiltered() && gt.isHet()){ //only phased Het SNP is considered
						GenomeLoc genomeLoc = glp.createGenomeLoc(vc.getContig(),vc.getStart(), vc.getEnd());
						Byte alleleA = gt.getAllele(0).getBases()[0];
						Byte alleleB = gt.getAllele(1).getBases()[0];
						Pair<Byte, Byte> tmp = new Pair<Byte, Byte>(alleleA, alleleB);
						genotypeMapList2.put(genomeLoc, tmp);
					}
				
				
			}
			itSnp.close();
			
			//read all reads within the region
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();	
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				
				if(failFlagFilter(r)){
					continue;
				}else{
					
					if(countedReadsEnd1.containsKey(r.getReadName())){
						countedReadsEnd1.remove(r.getReadName()); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd1.put(r.getReadName(), r);
					}
				}
				
			}
			wgsIt.close();
			
			HashMap<String, SAMRecord> countedReadsEnd2 = new HashMap<String, SAMRecord>();
			wgsIt = wgsReader.queryOverlapping(chr,start2,end2);
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				if(failFlagFilter(r)){
					continue;
				}else{
					if(countedReadsEnd2.containsKey(r.getReadName())){
						countedReadsEnd2.remove(r.getReadName()); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd2.put(r.getReadName(), r);
					}
				}
			}
			wgsIt.close();
			
			processBam(genotypeMapList1, genotypeMapList2, countedReadsEnd1, countedReadsEnd2, wgsReader, regionName, ob, exp, fdr);
		}
		
		private void processBam(HashMap<GenomeLoc, Pair<Byte, Byte>> genotypeMapList1, HashMap<GenomeLoc, Pair<Byte, Byte>> genotypeMapList2, HashMap<String, SAMRecord> countedReadsEnd1, 
				HashMap<String, SAMRecord> countedReadsEnd2, SamReader wgsReader, String regionName, double ob, double exp, double fdr){
			//end1
			for(GenomeLoc loc : genotypeMapList1.keySet()){
				Pair<Byte, Byte> alleles = genotypeMapList1.get(loc);
				SAMRecordIterator samIt = wgsReader.queryOverlapping(loc.getContig(), loc.getStart(), loc.getStop());

				Pair<Integer, Integer> ctAlleleRefAll = new Pair<Integer, Integer>(0,0);
				Pair<Integer, Integer> ctAlleleAltAll = new Pair<Integer, Integer>(0,0);
				Pair<Integer, Integer> ctAlleleRefEnd1 = new Pair<Integer, Integer>(0,0);
				Pair<Integer, Integer> ctAlleleAltEnd1 = new Pair<Integer, Integer>(0,0);
				Pair<Integer, Integer> ctAlleleRefEnd2 = new Pair<Integer, Integer>(0,0);
				Pair<Integer, Integer> ctAlleleAltEnd2 = new Pair<Integer, Integer>(0,0);
				
				int refAlleleCount = 0;
				int altAlleleCount = 0;


				while(samIt.hasNext()){
					SAMRecord r1 = samIt.next();
					boolean negStrand = r1.getReadNegativeStrandFlag();
					boolean secondPair = r1.getReadPairedFlag() && r1.getSecondOfPairFlag();

					
					byte[] readBases = r1.getReadBases();
					byte[] readBasesQ = r1.getBaseQualities();
					
					
					//define which allele the reads belonged to
					
					//in C/T or A/G SNP only consider bisulfite uncoverted strand
					if((SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.C) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.T))
							|| (SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.T) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.C))
							){
						if((!secondPair && !negStrand) || (secondPair && negStrand)){
							continue;
						}
					}else if((SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.G) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.A))
							|| (SequenceUtil.basesEqual(alleles.getFirst(), SequenceUtil.A) && SequenceUtil.basesEqual(alleles.getSecond(), SequenceUtil.G))){
						if((!secondPair && negStrand) || (secondPair && !negStrand)){
							continue;
						}
					}
					
					//int pos = loc.getStart() - d.getAlignmentStart();
					int pos = r1.getReadPositionAtReferencePosition(loc.getStart())-1;
					if(pos < 0) //no such a postion in the clipped reads
						continue;
						
					byte snpBase = readBases[pos];
					if(readBasesQ[pos] < minBaseQ){
						continue;
					}
					
					
					
					byte[] refBaseEnd1 = BisulfiteHicUtils.modifyRefSeqByCigar(refParser.loadFragment(r1.getAlignmentStart()-1, r1.getAlignmentEnd()-r1.getAlignmentStart()+1).getBytes(), r1.getCigarString());
					if(r1.getReadNegativeStrandFlag()){
						refBaseEnd1=BisulfiteHicUtils.complementArray(refBaseEnd1);
					}
					byte[] readBases1 = BisulfiteHicUtils.getClippedReadsBase(r1);
					//System.err.println(r1.getCigarString());
					int pos1 = 0;
					if(!turnOffBisulfiteFilter){
						pos1 = CcInferenceUtils.bisulfiteIncompleteReads(readBases1, refBaseEnd1, r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag(), false, true, true, 0.1, 0.4, 0, "CH", 1);
						if(pos1<0){
							continue;
						}
					}
					
					//System.err.println(r2.getCigarString());
					byte[] refBaseEnd2 = BisulfiteHicUtils.modifyRefSeqByCigar(refParser.loadFragment(r2.getAlignmentStart()-1, r2.getAlignmentEnd()-r2.getAlignmentStart()+1).getBytes(), r2.getCigarString());
					if(r2.getReadNegativeStrandFlag()){
						refBaseEnd2=BisulfiteHicUtils.complementArray(refBaseEnd2);
					}
					byte[] readBases2 = BisulfiteHicUtils.getClippedReadsBase(r2);
					int pos2 = 0;
					if(!turnOffBisulfiteFilter){
						pos2 = CcInferenceUtils.bisulfiteIncompleteReads(readBases2, refBaseEnd2, r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag(), false, true, false, 0.1, 0.4, 0, "CH", 1);
						if(pos2<0){
							continue;
						}
					}
					
					
					
					Pair<Integer, Integer> ctSummary = BisulfiteHicUtils.readsMethySummary(r1, minMapQ, minBaseQ);
					int numC = ctSummary.getFirst();
					int numT = ctSummary.getSecond();
					if((numC + numT) == 0)
						continue;
					if(BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getFirst(), snpBase, negStrand, secondPair)){
						int numCRef = ctAlleleRef.getFirst() + numC;
						int numTRef = ctAlleleRef.getSecond() + numT;
						ctAlleleRef = new Pair<Integer, Integer>(numCRef,numTRef);
						refAlleleCount++;
					}else if(BisulfiteHicUtils.baseMatchInBisulfiteSpace(alleles.getSecond(), snpBase, negStrand, secondPair)){
						int numCAlt = ctAlleleAlt.getFirst() + numC;
						int numTAlt = ctAlleleAlt.getSecond() + numT;
						ctAlleleAlt = new Pair<Integer, Integer>(numCAlt,numTAlt);
						altAlleleCount++;
					}
					
					
					
				}
				samIt.close();
				if(refAlleleCount < coverageRef || altAlleleCount< coverageAlt || (ctAlleleRefAll.getFirst() + ctAlleleRefAll.getSecond() ==0 && ctAlleleAltAll.getFirst() + ctAlleleAltAll.getSecond() ==0))
					continue;
				writer.println(loc.getContig() + "\t" + regionStart + "\t" + regionEnd + "\t" + loc.getStop()  + "\t" + String.format("%c", alleles.getFirst()) + "\t" + refAlleleCount
						  + "\t" + ctAlleleRef.getFirst() + "\t" + ctAlleleRef.getSecond() + "\t" +(double)ctAlleleRef.getFirst()/(double)(ctAlleleRef.getFirst() + ctAlleleRef.getSecond())
						  + "\t" + String.format("%c", alleles.getSecond()) + "\t" + altAlleleCount
						  + "\t" + ctAlleleAlt.getFirst() + "\t" + ctAlleleAlt.getSecond() + "\t" +(double)ctAlleleAlt.getFirst()/(double)(ctAlleleAlt.getFirst() + ctAlleleAlt.getSecond()));
				if(VCFCOUNT % PURGE_INTERVAL == 0 ){
					writer.flush();
					log.info("Processing allele ... " + VCFCOUNT);
				}
					
				VCFCOUNT++;
			}
			
		}
		
		private void methyCorrelation(String chr, int start1, int end1, int start2, int end2, SamReader wgsReader, double ob, double exp, double fdr, String bamFile) throws Exception{
			
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + start2 + ":" + end2;
			
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();	
			boolean reopenFlag=false;
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;
				if(queryNum % 1000000 == 0){
					reopenFlag=true;
				}
			//	if(r.getReadName().equalsIgnoreCase("K00168:28:H3NKLBBXX:7:1214:10064:40385")){
					//refParser.setCurrentSequence(chr);
					//byte[] refBaseEnd1 = refParser.loadFragment(r.getAlignmentStart()-1, r.getAlignmentEnd()-r.getAlignmentStart()+1).getBytes();
					//if(r.getReadNegativeStrandFlag()){
					//	refBaseEnd1=BisulfiteHicUtils.complementArray(refBaseEnd1);
					//}
					//byte[] readBases = BisulfiteHicUtils.getClippedReadsBase(r);
					//int pos = CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBaseEnd1, r.getReadNegativeStrandFlag(), r.getReadPairedFlag() & r.getSecondOfPairFlag(), false, true, false, 0.1, 0.4, 0, "CH", 1);
					//int truncateStart = r.getReadPositionAtReferencePosition(r.getAlignmentStart());
					//int truncateEnd = r.getReadPositionAtReferencePosition(r.getAlignmentEnd());
					
					//System.err.println(truncateStart + "\t" + truncateEnd );
					//System.err.println(new String(refBaseEnd1));
					//System.err.println(new String(readBases));
					//System.err.println(r.getReadName() + "\t" + countedReadsEnd1.containsKey(r.getReadName()) + "\t" + r + "\t" + pos);
				//}
				if(failFlagFilter(r)){
					continue;
				}else{
					//int bisulfiteStartPos = CcInferenceUtils.bisulfiteIncompleteReads(r);
					//if(bisulfiteStartPos<0){
					//	continue;
					//}
					//if(r.getReadName().equalsIgnoreCase("K00168:28:H3NKLBBXX:7:1214:10064:40385")){
					//	refParser.setCurrentSequence(chr);
					//	byte[] refBaseEnd1 = refParser.loadFragment(r.getAlignmentStart()-1, r.getAlignmentEnd()-r.getAlignmentStart()+1).getBytes();
					//	int pos = CcInferenceUtils.bisulfiteIncompleteReads(r.getReadBases(), refBaseEnd1, r.getReadNegativeStrandFlag(), r.getReadPairedFlag() & r.getSecondOfPairFlag(), false, true, false, 0.1, 0.4, 0, "CH", 1);
						
					//	System.err.println(new String(refBaseEnd1));
					//	System.err.println(new String(r.getReadBases()));
					//	System.err.println(r.getReadName() + "\t" + countedReadsEnd1.containsKey(r.getReadName()) + "\t" + r + "\t" + pos);
					//}
					if(countedReadsEnd1.containsKey(r.getReadName())){
						countedReadsEnd1.remove(r.getReadName()); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd1.put(r.getReadName(), r);
					}
				}
				
			}
			wgsIt.close();
			if(reopenFlag){
				//wgsReader.close();
				//wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
				reopenFlag=false;
			}
			HashMap<String, SAMRecord> countedReadsEnd2 = new HashMap<String, SAMRecord>();
			wgsIt = wgsReader.queryOverlapping(chr,start2,end2);
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;
				if(queryNum % 1000000 == 0){
					reopenFlag=true;
				}
				//if(r.getReadName().equalsIgnoreCase("K00168:28:H3NKLBBXX:7:1214:10064:40385")){
					//byte[] refBaseEnd1 = refParser.loadFragment(r.getAlignmentStart()-1, r.getAlignmentEnd()-r.getAlignmentStart()+1).getBytes();
					//if(r.getReadNegativeStrandFlag()){
					//	refBaseEnd1=BisulfiteHicUtils.complementArray(refBaseEnd1);
					//}
					//byte[] readBases = BisulfiteHicUtils.getClippedReadsBase(r);
					//int truncateStart = r.getReadPositionAtReferencePosition(r.getAlignmentStart());
					//int truncateEnd = r.getReadPositionAtReferencePosition(r.getAlignmentEnd());
					
					//System.err.println(truncateStart + "\t" + truncateEnd );
					//System.err.println(new String(refBaseEnd1));
					//System.err.println(new String(readBases));
					
					//int pos = CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBaseEnd1, r.getReadNegativeStrandFlag(), r.getReadPairedFlag() & r.getSecondOfPairFlag(), false, true, false, 0.1, 0.4, 0, "CH", 1);
					//System.err.println(r.getReadName() + "\t" + countedReadsEnd1.containsKey(r.getReadName()) + "\t" + r + "\t" + pos);
					//System.err.println(new String(refBaseEnd1));
					//System.err.println(new String(readBases));
				//}
				if(failFlagFilter(r)){
					continue;
				}else{
					//int bisulfiteStartPos = CcInferenceUtils.bisulfiteIncompleteReads(r);
					//if(bisulfiteStartPos<0){
					//	continue;
					//}
					//if(r.getReadName().equalsIgnoreCase("K00168:28:H3NKLBBXX:7:1214:10064:40385")){
					//	System.err.println(r.getReadName() + "\t" + countedReadsEnd2.containsKey(r.getReadName()) + "\t" + r);
					//}
					if(countedReadsEnd2.containsKey(r.getReadName())){
						countedReadsEnd2.remove(r.getReadName()); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd2.put(r.getReadName(), r);
					}
				}
			}
			wgsIt.close();
			if(reopenFlag){
				//wgsReader.close();
				//wgsReader = SamReaderFactory.makeDefault().setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true).validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
				log.info("reopen bam file to avoid too many connections...");
				reopenFlag=false;
			}
			//System.err.println(countedReadsEnd1.size() + "\t" + countedReadsEnd2.size());
			int readNum = 0;
			for(String readName : countedReadsEnd1.keySet()){
				if(countedReadsEnd2.containsKey(readName)){
					
					SAMRecord r1 = countedReadsEnd1.get(readName);
					SAMRecord r2 = countedReadsEnd2.get(readName);
					
					if(!chr.equalsIgnoreCase(preChr)){
						refParser.close();
						refParser.setCurrentSequence(chr);
						preChr = chr;
					}
					
					byte[] refBaseEnd1 = BisulfiteHicUtils.modifyRefSeqByCigar(refParser.loadFragment(r1.getAlignmentStart()-1, r1.getAlignmentEnd()-r1.getAlignmentStart()+1).getBytes(), r1.getCigarString());
					if(r1.getReadNegativeStrandFlag()){
						refBaseEnd1=BisulfiteHicUtils.complementArray(refBaseEnd1);
					}
					byte[] readBases1 = BisulfiteHicUtils.getClippedReadsBase(r1);
					//System.err.println(r1.getCigarString());
					int pos1 = 0;
					if(!turnOffBisulfiteFilter){
						pos1 = CcInferenceUtils.bisulfiteIncompleteReads(readBases1, refBaseEnd1, r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag(), false, true, true, 0.1, 0.4, 0, "CH", 1);
						if(pos1<0){
							continue;
						}
					}
					
					//System.err.println(r2.getCigarString());
					byte[] refBaseEnd2 = BisulfiteHicUtils.modifyRefSeqByCigar(refParser.loadFragment(r2.getAlignmentStart()-1, r2.getAlignmentEnd()-r2.getAlignmentStart()+1).getBytes(), r2.getCigarString());
					if(r2.getReadNegativeStrandFlag()){
						refBaseEnd2=BisulfiteHicUtils.complementArray(refBaseEnd2);
					}
					byte[] readBases2 = BisulfiteHicUtils.getClippedReadsBase(r2);
					int pos2 = 0;
					if(!turnOffBisulfiteFilter){
						pos2 = CcInferenceUtils.bisulfiteIncompleteReads(readBases2, refBaseEnd2, r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag(), false, true, false, 0.1, 0.4, 0, "CH", 1);
						if(pos2<0){
							continue;
						}
					}
					
					Pair<Double,Integer> methySummaryEnd1 = BisulfiteHicUtils.readsMethySummaryWithFivePrimeFilter(refBaseEnd1,readBases1, BisulfiteHicUtils.getClippedReadsBaseQuality(r1),
														r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag(), minBaseQ, pos1);
					Pair<Double,Integer> methySummaryEnd2 = BisulfiteHicUtils.readsMethySummaryWithFivePrimeFilter(refBaseEnd2,readBases2, BisulfiteHicUtils.getClippedReadsBaseQuality(r2),
							r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag(), minBaseQ, pos2);
					//chr1,read.start1,read.end1,chr2,read.start2,read.end2,region_name, id, hiccups_ob, hiccups_exp, hiccups_fdr, methy_end1, cpg_num_end1, methy_end2, cpg_num_end2
					if(!Double.isNaN(methySummaryEnd1.getFirst()) && !Double.isNaN(methySummaryEnd2.getFirst())){
						methyCorWriter.write(chr + "\t" + start1 + "\t" + end1 + "\t" + chr + "\t" + start2 + "\t" + end2 + "\t" + 
								regionName + "\t" + readNum + "\t" + ob+ "\t" + exp+ "\t" + fdr + "\t" + methySummaryEnd1.getFirst() + "\t" + methySummaryEnd1.getSecond()
								 + "\t" + methySummaryEnd2.getFirst() + "\t" + methySummaryEnd2.getSecond() + "\n");
								readNum++;
					}
					
				}
				
				
			}
			//System.err.println(readNum);
		}
		
		private void methyCorrelationRandom(String chr, int start1, int end1, int start2, int end2, SamReader wgsReader, double ob, double exp, double fdr, String bamFile) throws Exception{
			
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + start2 + ":" + end2;
			
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();	
			boolean reopenFlag=false;
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;
				if(queryNum % 1000000 == 0){
					reopenFlag=true;
				}
			
				if(failFlagFilter(r)){
					continue;
				}else{
					
					if(countedReadsEnd1.containsKey(r.getReadName())){
						countedReadsEnd1.remove(r.getReadName()); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd1.put(r.getReadName(), r);
					}
				}
				
			}
			wgsIt.close();
			
			HashMap<String, SAMRecord> countedReadsEnd2 = new HashMap<String, SAMRecord>();
			wgsIt = wgsReader.queryOverlapping(chr,start2,end2);
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				queryNum++;
				if(queryNum % 1000000 == 0){
					reopenFlag=true;
				}
				
				if(failFlagFilter(r)){
					continue;
				}else{
					
					if(countedReadsEnd2.containsKey(r.getReadName())){
						countedReadsEnd2.remove(r.getReadName()); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd2.put(r.getReadName(), r);
					}
				}
			}
			wgsIt.close();
			
			//System.err.println(countedReadsEnd1.size() + "\t" + countedReadsEnd2.size());
			int readNum = 0;
			String[] readNames1 = countedReadsEnd1.keySet().toArray(new String[countedReadsEnd1.keySet().size()]);
			String[] readNames2 = countedReadsEnd2.keySet().toArray(new String[countedReadsEnd2.keySet().size()]);
			if(readNames1.length>10 && readNames2.length>10){
				for(String readName : readNames1){
					SAMRecord r1 = countedReadsEnd1.get(readName);
					int randomItem = generator.nextInt(readNames2.length);
					while(readNames2[randomItem].equalsIgnoreCase(readName)){
						randomItem = generator.nextInt(readNames2.length);
					}
					SAMRecord r2=countedReadsEnd2.get(readNames2[randomItem]);
					
					
					if(!chr.equalsIgnoreCase(preChr)){
						refParser.close();
						refParser.setCurrentSequence(chr);
						preChr = chr;
					}
					
					byte[] refBaseEnd1 = BisulfiteHicUtils.modifyRefSeqByCigar(refParser.loadFragment(r1.getAlignmentStart()-1, r1.getAlignmentEnd()-r1.getAlignmentStart()+1).getBytes(), r1.getCigarString());
					if(r1.getReadNegativeStrandFlag()){
						refBaseEnd1=BisulfiteHicUtils.complementArray(refBaseEnd1);
					}
					byte[] readBases1 = BisulfiteHicUtils.getClippedReadsBase(r1);
					//System.err.println(r1.getCigarString());
					int pos1 = 0;
					if(!turnOffBisulfiteFilter){
						pos1 = CcInferenceUtils.bisulfiteIncompleteReads(readBases1, refBaseEnd1, r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag(), false, true, true, 0.1, 0.4, 0, "CH", 1);
						if(pos1<0){
							continue;
						}
					}
					
					//System.err.println(r2.getCigarString());
					byte[] refBaseEnd2 = BisulfiteHicUtils.modifyRefSeqByCigar(refParser.loadFragment(r2.getAlignmentStart()-1, r2.getAlignmentEnd()-r2.getAlignmentStart()+1).getBytes(), r2.getCigarString());
					if(r2.getReadNegativeStrandFlag()){
						refBaseEnd2=BisulfiteHicUtils.complementArray(refBaseEnd2);
					}
					byte[] readBases2 = BisulfiteHicUtils.getClippedReadsBase(r2);
					int pos2 = 0;
					if(!turnOffBisulfiteFilter){
						pos2 = CcInferenceUtils.bisulfiteIncompleteReads(readBases2, refBaseEnd2, r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag(), false, true, false, 0.1, 0.4, 0, "CH", 1);
						if(pos2<0){
							continue;
						}
					}
					
					Pair<Double,Integer> methySummaryEnd1 = BisulfiteHicUtils.readsMethySummaryWithFivePrimeFilter(refBaseEnd1,readBases1, BisulfiteHicUtils.getClippedReadsBaseQuality(r1),
														r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag(), minBaseQ, pos1);
					Pair<Double,Integer> methySummaryEnd2 = BisulfiteHicUtils.readsMethySummaryWithFivePrimeFilter(refBaseEnd2,readBases2, BisulfiteHicUtils.getClippedReadsBaseQuality(r2),
							r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag(), minBaseQ, pos2);
					//chr1,read.start1,read.end1,chr2,read.start2,read.end2,region_name, id, hiccups_ob, hiccups_exp, hiccups_fdr, methy_end1, cpg_num_end1, methy_end2, cpg_num_end2
					if(!Double.isNaN(methySummaryEnd1.getFirst()) && !Double.isNaN(methySummaryEnd2.getFirst())){
						methyCorWriter.write(chr + "\t" + start1 + "\t" + end1 + "\t" + chr + "\t" + start2 + "\t" + end2 + "\t" + 
								regionName + "\t" + readNum + "\t" + ob+ "\t" + exp+ "\t" + fdr + "\t" + methySummaryEnd1.getFirst() + "\t" + methySummaryEnd1.getSecond()
								 + "\t" + methySummaryEnd2.getFirst() + "\t" + methySummaryEnd2.getSecond() + "\n");
								readNum++;
					}
					if(randomSampleReads>0 && readNum > randomSampleReads){
						break;
					}
					
					
				}
			}
			

			//System.err.println(readNum);
		}
		
		private void methyCorrelationShort(String chr, int start1, int end1, int start2, int end2, SamReader wgsReader, double ob, double exp, double fdr, String bamFile) throws Exception{
			
			String regionName = chr + ":" + start1 + ":" + end1 + ":" + start2 + ":" + end2;
			
			SAMRecordIterator wgsIt = wgsReader.queryOverlapping(chr,start1,end1);
			HashMap<String, SAMRecord> countedReadsEnd1 = new HashMap<String, SAMRecord>();	
			HashMap<String, SAMRecord> countedReadsEnd2 = new HashMap<String, SAMRecord>();
			

			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				
				if(failFlagFilter(r) || Math.abs(r.getInferredInsertSize())>onlyShort){
					continue;
				}else{
					
					if(countedReadsEnd1.containsKey(r.getReadName())){
						countedReadsEnd2.put(r.getReadName(), r); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd1.put(r.getReadName(), r);
					}
				}
				
			}
			wgsIt.close();

			wgsIt = wgsReader.queryOverlapping(chr,start2,end2);
			while(wgsIt.hasNext()){
				SAMRecord r = wgsIt.next();
				if(failFlagFilter(r) || Math.abs(r.getInferredInsertSize())>onlyShort){
					continue;
				}else{
					if(countedReadsEnd1.containsKey(r.getReadName())){
						countedReadsEnd2.put(r.getReadName(), r); //if two ends are within the same half end of hiccup links
					}else{
						countedReadsEnd1.put(r.getReadName(), r);
					}
				}
			}
			wgsIt.close();
			
			//System.err.println(countedReadsEnd1.size() + "\t" + countedReadsEnd2.size());
			int readNum = 0;
			for(String readName : countedReadsEnd1.keySet()){
				if(countedReadsEnd2.containsKey(readName)){
					
					SAMRecord r1 = countedReadsEnd1.get(readName);
					SAMRecord r2 = countedReadsEnd2.get(readName);
					
					if(!chr.equalsIgnoreCase(preChr)){
						refParser.close();
						refParser.setCurrentSequence(chr);
						preChr = chr;
					}
					
					byte[] refBaseEnd1 = BisulfiteHicUtils.modifyRefSeqByCigar(refParser.loadFragment(r1.getAlignmentStart()-1, r1.getAlignmentEnd()-r1.getAlignmentStart()+1).getBytes(), r1.getCigarString());
					if(r1.getReadNegativeStrandFlag()){
						refBaseEnd1=BisulfiteHicUtils.complementArray(refBaseEnd1);
					}
					byte[] readBases1 = BisulfiteHicUtils.getClippedReadsBase(r1);
					//System.err.println(r1.getCigarString());
					int pos1 = 0;
					if(!turnOffBisulfiteFilter){
						pos1 = CcInferenceUtils.bisulfiteIncompleteReads(readBases1, refBaseEnd1, r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag(), false, true, true, 0.1, 0.4, 0, "CH", 1);
						if(pos1<0){
							continue;
						}
					}
					
					//System.err.println(r2.getCigarString());
					byte[] refBaseEnd2 = BisulfiteHicUtils.modifyRefSeqByCigar(refParser.loadFragment(r2.getAlignmentStart()-1, r2.getAlignmentEnd()-r2.getAlignmentStart()+1).getBytes(), r2.getCigarString());
					if(r2.getReadNegativeStrandFlag()){
						refBaseEnd2=BisulfiteHicUtils.complementArray(refBaseEnd2);
					}
					byte[] readBases2 = BisulfiteHicUtils.getClippedReadsBase(r2);
					int pos2 = 0;
					if(!turnOffBisulfiteFilter){
						pos2 = CcInferenceUtils.bisulfiteIncompleteReads(readBases2, refBaseEnd2, r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag(), false, true, false, 0.1, 0.4, 0, "CH", 1);
						if(pos2<0){
							continue;
						}
					}
					
					Pair<Double,Integer> methySummaryEnd1 = BisulfiteHicUtils.readsMethySummaryWithFivePrimeFilter(refBaseEnd1,readBases1, BisulfiteHicUtils.getClippedReadsBaseQuality(r1),
														r1.getReadNegativeStrandFlag(), r1.getReadPairedFlag() & r1.getSecondOfPairFlag(), minBaseQ, pos1);
					Pair<Double,Integer> methySummaryEnd2 = BisulfiteHicUtils.readsMethySummaryWithFivePrimeFilter(refBaseEnd2,readBases2, BisulfiteHicUtils.getClippedReadsBaseQuality(r2),
							r2.getReadNegativeStrandFlag(), r2.getReadPairedFlag() & r2.getSecondOfPairFlag(), minBaseQ, pos2);
					//chr1,read.start1,read.end1,chr2,read.start2,read.end2,region_name, id, hiccups_ob, hiccups_exp, hiccups_fdr, methy_end1, cpg_num_end1, methy_end2, cpg_num_end2
					if(!Double.isNaN(methySummaryEnd1.getFirst()) && !Double.isNaN(methySummaryEnd2.getFirst())){
						methyCorWriter.write(chr + "\t" + start1 + "\t" + end1 + "\t" + chr + "\t" + start2 + "\t" + end2 + "\t" + 
								regionName + "\t" + readNum + "\t" + ob+ "\t" + exp+ "\t" + fdr + "\t" + methySummaryEnd1.getFirst() + "\t" + methySummaryEnd1.getSecond()
								 + "\t" + methySummaryEnd2.getFirst() + "\t" + methySummaryEnd2.getSecond() + "\n");
								readNum++;
					}
					if(randomSampleReads>0 && readNum > randomSampleReads){
						break;
					}
					
				}
				
				
			}
			//System.err.println(readNum);
		}
		
		private boolean failFlagFilter(SAMRecord r){
			return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
					|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || !r.getProperPairFlag();
		}

		
		
		
		private void initiate(String outputPrefix, String refFile) throws Exception{
			startTime = System.currentTimeMillis();
			refParser = new TwoBitParser(new File(refFile));
			
			if(vcf != null){
				asmWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".asm.txt.gz")), "UTF-8");
				File indexFile = Tribble.tabixIndexFile(new File(vcf));
				if(!indexFile.exists() || !indexFile.canRead()){
					throw new UnsortedFileException(vcf + " file's index " + indexFile.getName() + " does not exist, please use tabix to index it ...");
				}

				vcfReader = new VCFFileReader(new File(vcf), true);
			}else{
				methyCorWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".methyCor.txt.gz")), "UTF-8");
			}
			
			generator = new MersenneTwister();
		}

		private void finish() throws Exception{
			refParser.closeParser();
			
			
			if(vcf != null){
				vcfReader.close();
				asmWriter.close();
			}else{
				methyCorWriter.close();
			}
			
			
			long endTime   = System.currentTimeMillis();
			double totalTime = endTime - startTime;
			totalTime /= 1000;
			double totalTimeMins = totalTime/60;
			double totalTimeHours = totalTime/3600;
			
			log.info("MethyCorAcrossHiccups's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
		}
		

}
