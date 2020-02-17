# Bhmem
A mapping tool for Methyl-HiC and single-cell Methyl-HiC based on bwa


## pre-requisite:
* JAVA oracle JDK8
* GNU Make 3.81
* gcc 4.8.2
* wget
* maven

## command to install it
```
unzip dnaase-bisulfitehic-5dbb6ce7bc3c.zip
cd dnaase-bisulfitehic-5dbb6ce7bc3c
source ./install.sh
conda create -n bisulfitehic python=3.6
source activate bisulfitehic
pip install numpy

export JAVA_HOME=/usr/libexec/java_home
```


## prepare index file for genome and bisulfite converted genomes
* Use bismark_genome_preparation in bismark to prepare the bisulfite converted genomes. It will generate Bisulfite_Genome directory under the same directory as your reference genome, such as mm9.fa.

* Then use bwa index to generate index for each of the converted genomes (Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa and Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa)

* under your reference genome folder, you should have mm9.fa, mm9.fa.fai, mm9.dict, Bisulfite_Genome (directory), Bisulfite_Genome/CT_conversion (directory),  Bisulfite_Genome/GA_conversion (directory), Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa (and its bwa indexed files), and Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa (and its bwa indexed files)

## test
```
./bisulfitehicMap /Users/yaping/Documents/workspace/data/genomic_data/genome/mm9/mm9.fa ~/Downloads/test.bam /Users/yaping/Documents/workspace/project/bisulfiteHaploseq/testData/bisulfite/test_R1.fq.gz /Users/yaping/Documents/workspace/project/bisulfiteHaploseq/testData/bisulfite/test_R2.fq.gz  -rgId test -rgSm test1
```

## test on standalone jar file:
```
java -Xmx5G -Djava.library.path=/Users/yaping/Documents/workspace/code/jbwa/src/main/native/ -cp "target/bisulfitehic-0.01-jar-with-dependencies.jar:lib/jbwa.jar" main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem /Users/yaping/Documents/workspace/data/genomic_data/genome/mm9/mm9.fa /Users/yaping/Documents/workspace/project/bisulfiteHaploseq/testData/bisulfite/testOut.bam /Users/yaping/Documents/workspace/project/bisulfiteHaploseq/testData/bisulfite/test_R1.fq.gz /Users/yaping/Documents/workspace/project/bisulfiteHaploseq/testData/bisulfite/test_R2.fq.gz
```

```
java -Xmx5G -jar target/bisulfitehic-0.01-jar-with-dependencies.jar -Djava.library.path=/Users/yaping/Documents/workspace/code/jbwa/src/main/native/ -cp "lib/jbwa.jar" /Users/yaping/Documents/workspace/data/genomic_data/genome/mm9/mm9.fa /Users/yaping/Documents/workspace/project/bisulfiteHaploseq/testData/bisulfite/testOut.bam /Users/yaping/Documents/workspace/project/bisulfiteHaploseq/testData/bisulfite/test_R1.fq.gz /Users/yaping/Documents/workspace/project/bisulfiteHaploseq/testData/bisulfite/test_R2.fq.gz
```










