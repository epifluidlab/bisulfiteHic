/**
 * enhanceHicMatrix.java
 * Nov 8, 2016
 * 9:26:51 AM
 * yaping    lyping1986@gmail.com
 */
package main.java.edu.mit.compbio.bisulfitehic.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import main.java.edu.mit.compbio.ccinference.utils.CcInferenceUtils;

import org.apache.commons.collections4.IteratorUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * use the concordance of the methylation status at each end to increase the resolution of HiC. 
 * output as bed format (gzipped):
 * chr6    80000_120000    120000_160000   10.278
 *chr6    80000_120000    160000_200000   3.648
 *chr6    80000_120000    200000_240000   4.204
 *
 * <readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>
 * 
 */
public class SingleBaseResHicMatrix {

	
	@Option(name="-minBaseQ",usage="minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;
	
	@Option(name="-turnOffBisulfiteFilter",usage="turn off bisulfite conversion reads filter. Default: false")
	public boolean turnOffBisulfiteFilter = false;
	
	@Option(name="-minMapQ",usage="minimum mapping quality score required to check. Default: 30")
	public int minMapQ = -1;

	@Option(name="-minDistEnds",usage="minimum distance between two end to be output. Default: 20000")
	public int minDistEnds = 20000;

	@Option(name="-maxDistEnds",usage="maximum distance between two end to be output. Default: 3000000")
	public int maxDistEnds = 3000000;

	
	//@Option(name="-minCgNum",usage="For ASM mode, specify the minimum  number of CpG to call ASM. Default: 1")
	//public int minCgNum = 1;
	
	@Option(name="-permutations",usage="number of permutations to generate expectated p value within the same segments, default: 10")
	public int permutations = 10;
	
	@Option(name="-skipProcessFirstRow",usage="don't process first row, default: not enabled")
	public boolean skipProcessFirstRow = false;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "SingleBaseResHicMatrix [opts] outputPrefix cpg_list.bed[.gz] input.bam";
	
	private static Logger log = Logger.getLogger(SingleBaseResHicMatrix.class);

	private OutputStreamWriter writer = null; 
	private MersenneTwister generator = null;
	
	
	private static long startTime = -1;
	private static long queryNum=0;

	/**
	 * @param args
	 */
	public static void main(String[] args)
		throws Exception {
			SingleBaseResHicMatrix ehm = new SingleBaseResHicMatrix();
				BasicConfigurator.configure();
				ehm.doMain(args);
	}
			

	public void doMain(String[] args)
	throws Exception {

		CmdLineParser parser = new CmdLineParser(this);
							//parser.setUsageWidth(80);
							try
							{
								if(help || args.length < 3) throw new CmdLineException(parser, USAGE, new Throwable());
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
							String outputPrefix = arguments.get(0);
							String cpgFile = arguments.get(1);
							String bamFile = arguments.get(2);

							initiate(outputPrefix);
							
							log.info("Parsing input cpgList bed file ...");
							
							
							HashMap<String,IntervalTree<String>> allCpgLocCollections = new HashMap<String,IntervalTree<String>>();
								log.info("Loading all CpG intervals ... ");
								GZIPInputStream gzipInputStream = null;
								BufferedReader br1;
								if(cpgFile.endsWith(".gz")){
									gzipInputStream = new GZIPInputStream(new FileInputStream(cpgFile));
									br1 = new BufferedReader(new InputStreamReader(gzipInputStream));
									
								}else{
									br1 = new BufferedReader(new FileReader(cpgFile));
								}
									
									String line1;
									
									while( (line1 = br1.readLine()) != null){
										if(line1.startsWith("#"))
											continue;
										String[] splitLines = line1.split("\t");
										if(splitLines.length<3){
											continue;
										}
										String chr = splitLines[0];
										int start = Integer.parseInt(splitLines[1]);
										int end = Integer.parseInt(splitLines[2]);
										IntervalTree<String> tree;
										
										if(allCpgLocCollections.containsKey(chr)){
											tree = allCpgLocCollections.get(chr);
										}else{
											tree = new IntervalTree<String>();
										}
										String strand = ".";
										if(splitLines.length >= 6){
											if(splitLines[5].equalsIgnoreCase("-")){
												strand = "-";
											}else if(splitLines[5].equalsIgnoreCase("+")){
												strand = "+";
											}
											//strand = splitLines[5].equalsIgnoreCase("-") ? "-" : "+";
										}
										tree.put(start, end, strand);
										allCpgLocCollections.put(chr, tree);
										
										
									}
									if(cpgFile.endsWith(".gz")){
										gzipInputStream.close();
									}
									br1.close();
							
					log.info("Extracting reads within intervals that contains cpg sites ... ");
					HashMap<String,IntervalTree<HashMap<String,Integer>>> cpgMethyMatrix = new HashMap<String,IntervalTree<HashMap<String,Integer>>>();
					
					SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
					
					SAMRecordIterator wgsIt = wgsReader.iterator();
					
					while(wgsIt.hasNext()){
						SAMRecord r = wgsIt.next();
						String readName = r.getReadName();
						String chr = r.getContig();
						int start = r.getAlignmentStart();
						int end = r.getAlignmentEnd();
						Iterator<IntervalTree.Node<String>> cpgs = null;
						if(allCpgLocCollections.containsKey(chr)){
							cpgs = allCpgLocCollections.get(chr).overlappers(start, end);
							if(cpgs == null){
								continue;
							}
							//if(allCpgLocCollections.get(chr).minOverlapper(start, end)==null){
							//	continue;
							//}else{
							//	
							//}
							
						}else{
							continue;
						}
						
						if(failFlagFilter(r)){
							continue;
						}else{
							TreeMap<Integer, Integer> cpgSummary = BisulfiteHicUtils.readsMethySummaryAtEachLoc(r, minMapQ, minBaseQ, turnOffBisulfiteFilter);
							while(cpgs.hasNext()){
								IntervalTree.Node<String> cpg = cpgs.next();
								int pos = cpg.getStart() - start;
								if(cpgSummary.containsKey(pos)){
									IntervalTree<HashMap<String,Integer>> tmpInterval = null;
									if(cpgMethyMatrix.containsKey(chr)){
										tmpInterval = cpgMethyMatrix.get(chr);
										HashMap<String,Integer> tmpRead = tmpInterval.find(cpg.getStart(), cpg.getEnd()).getValue();
										if(tmpRead == null){
											tmpRead = new HashMap<String,Integer>();
											tmpRead.put(readName, cpgSummary.get(pos));
										}else{
											if(tmpRead.containsKey(readName)){
												// two ends overlapped and have not be trimmed?
											}else{
												tmpRead.put(readName, cpgSummary.get(pos));
											}
										}
										tmpInterval.put(cpg.getStart(), cpg.getEnd(), tmpRead);
									}else{
										tmpInterval = new IntervalTree<HashMap<String,Integer>>();
										HashMap<String,Integer> tmpRead = new HashMap<String,Integer>();
										tmpRead.put(readName, cpgSummary.get(pos));
										tmpInterval.put(cpg.getStart(), cpg.getEnd(), tmpRead);
									}
									cpgMethyMatrix.put(chr, tmpInterval);
									
								}
							}
						
						}

					}
					wgsIt.close();
					wgsReader.close();
					
					log.info("Processing each pair of CpGs' interactin frequency ... ");
					for(String chr : cpgMethyMatrix.keySet()){
						IntervalTree<HashMap<String,Integer>> cpgs = cpgMethyMatrix.get(chr);
						List<Node<HashMap<String, Integer>>> cpgList  = IteratorUtils.toList(cpgs.iterator());  
						for(int i = 0; i < cpgList.size()-1; i++){
							Node<HashMap<String, Integer>> cpgi = cpgList.get(i);
							for(int j = i+1; j < cpgList.size(); j++){
								Node<HashMap<String, Integer>> cpgj = cpgList.get(j);
								if(Math.abs(cpgi.getEnd() - cpgj.getEnd()) > maxDistEnds || Math.abs(cpgi.getEnd() - cpgj.getEnd()) < minDistEnds){
									continue;
								}

								ArrayList<Double> ratios = cpgContactFreq(cpgi.getValue(), cpgj.getValue());
								queryNum++;
								if(!ratios.get(0).isNaN()){
									//writer.write(chr + "\t" + cpgi.getStart() + "_" + cpgi.getEnd() + "\t" + cpgj.getStart() + "_" + cpgj.getEnd() + "\t" + ratios.getFirst() + "\t" + ratios.getSecond());
									String name = chr + ":" + cpgi.getStart() + ":" + cpgi.getEnd() + ":" + cpgj.getStart() + ":" + cpgj.getEnd();
									writer.write(chr + "\t" + cpgi.getStart() + "\t" + cpgi.getEnd() + "\t" + cpgj.getStart() + "\t" + cpgj.getEnd() + "\t" + name);
									for(Double ratio : ratios){
										writer.write("\t" + ratio);
									}
									writer.write("\n");
								}
								if(queryNum % 1000000 == 0){
									log.info("Processing CpG ... " + queryNum);
									writer.flush();
								}
							}
						}
						
					}
					
					
					finish();
			}
			
	//observed_ratio, permutated_ratio, correlation, p value, number of points, cor_random, pvalue_random, numOfPoints_random
			private ArrayList<Double> cpgContactFreq(HashMap<String, Integer> x, HashMap<String, Integer> y){
				int concordCount = 0;
				int disconcordCount = 0;
				ArrayList<Double> results = new ArrayList<Double>();
				
				ArrayList<Double> xList = new ArrayList<Double>();
				ArrayList<Double> yList = new ArrayList<Double>();
				
				for(String readName : x.keySet()){
					if(y.containsKey(readName)){
						if(x.get(readName) == y.get(readName)){
							concordCount++;
						}else{
							disconcordCount++;
						}
						xList.add((double)x.get(readName));
						yList.add((double)y.get(readName));
					}
				}
				
				if(concordCount + disconcordCount == 0){
					results.add(Double.NaN);
					results.add(Double.NaN);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					results.add(0.);
					
					results.add(Double.NaN);
					results.add(Double.NaN);
					results.add(0.);
					
					return results;
				}
				double observedRatio = Math.log10(concordCount==0 ? concordCount + 0.0001 : concordCount) - Math.log10(disconcordCount==0 ? disconcordCount + 0.0001 : disconcordCount);
				
				double[] xArray = ArrayUtils.toPrimitive(xList.toArray(new Double[xList.size()]));
				double[] yArray = ArrayUtils.toPrimitive(yList.toArray(new Double[yList.size()]));
				double[][] d = new double[xArray.length][2];
				for(int i = 0; i < xArray.length; i++){
					d[i][0] = xArray[i];
					d[i][1] = yArray[i];
				}
				PearsonsCorrelation corObserved = new PearsonsCorrelation(d);
				Double corEff = corObserved.correlation(xArray, yArray);
				Double pvalue = Double.NaN;
				if(!Double.isNaN(corEff)){
					pvalue = corObserved.getCorrelationPValues().getEntry(0, 1);
				}
				
				
				
				
				//permutation for the expectations:
				ArrayList<Double> xListRandom = new ArrayList<Double>();
				ArrayList<Double> yListRandom = new ArrayList<Double>();
				
				double expectRatio = 0.;
				List<String> keysY  = new ArrayList<String>(y.keySet());
				for(int i = 0; i < permutations; i++){
					int concordCountRandom = 0;
					int disconcordCountRandom = 0;
					
					for(String readName : x.keySet()){
						String randomKey = keysY.get( generator.nextInt(keysY.size()));
						if(x.get(readName) == y.get(randomKey)){
							concordCountRandom++;
						}else{
							disconcordCountRandom++;
						}
						xListRandom.add((double)x.get(readName));
						yListRandom.add((double)y.get(randomKey));
					}
					
					if(concordCountRandom + disconcordCountRandom == 0){
						continue;
					}
					expectRatio += Math.log10(concordCountRandom==0 ? concordCountRandom + 0.0001 : concordCountRandom) - Math.log10(disconcordCountRandom==0 ? disconcordCountRandom + 0.0001 : disconcordCountRandom);
				}
				
				expectRatio /= permutations;
				
				double[] xArrayRandom = ArrayUtils.toPrimitive(xListRandom.toArray(new Double[xList.size()]));
				double[] yArrayRandom = ArrayUtils.toPrimitive(yListRandom.toArray(new Double[yList.size()]));
				double[][] dRandom = new double[xArrayRandom.length][2];
				for(int j = 0; j < xArrayRandom.length; j++){
					dRandom[j][0] = xArrayRandom[j];
					dRandom[j][1] = yArrayRandom[j];
				}
				PearsonsCorrelation corObservedRandom = new PearsonsCorrelation(dRandom);
				Double corEffRandom = corObservedRandom.correlation(xArrayRandom, yArrayRandom);
				Double pvalueRandom = Double.NaN;
				if(!Double.isNaN(corEffRandom)){
					pvalueRandom = corObservedRandom.getCorrelationPValues().getEntry(0, 1);
				}

				results.add(observedRatio);
				results.add(expectRatio);
				
				results.add(corEff);
				results.add(pvalue);
				results.add((double)xArray.length);
				
				results.add(corEffRandom);
				results.add(pvalueRandom);
				results.add((double)xArrayRandom.length);
				
				
				return results;
			}
	
			
			private boolean failFlagFilter(SAMRecord r){
				return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
						|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag() || !r.getProperPairFlag() ;
			}

			
			
			
			private void initiate(String outputPrefix) throws Exception{
				startTime = System.currentTimeMillis();
				writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputPrefix + ".singleBaseResHicMatrix.bed.gz")), "UTF-8");
				
				generator = new MersenneTwister();
			}

			private void finish() throws Exception{

				writer.close();

				long endTime   = System.currentTimeMillis();
				double totalTime = endTime - startTime;
				totalTime /= 1000;
				double totalTimeMins = totalTime/60;
				double totalTimeHours = totalTime/3600;
				
				log.info("SingleBaseResHicMatrix's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
			}
			

}
