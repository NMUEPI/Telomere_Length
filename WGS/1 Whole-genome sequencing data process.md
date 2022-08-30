#Whole-genome sequencing data process

####*Align sequence reads to genome reference*

##1 Trimmomatic
```java
mkdir log_path

java -Xmx10g -jar trimmomatic-0.36.jar \
PE -phred33 -threads ${cpus} \
${sample}.R1.fastq.gz \
${sample}.R2.fastq.gz \
${sample}_R1_trim.fq.gz ${sample}_R1_trim_up.fq.gz \
${sample}_R2_trim.fq.gz ${sample}_R2_trim_up.fq.gz \
ILLUMINACLIP:${params.trim_adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 \
2> log_path/${sample}-trim.log
```


##2.1 BWA
```Shell
mkdir log_path

${bwa} mem -t ${cpus} \
	-R "@RG\tID:${sample}\tPL:illumina\tPU:${sample}_LCB_WGS\tSM:${sample}" \
	GRCh37.fa ${sample}_R1_trim.fq.gz ${sample}_R2_trim.fq.gz \
	| ${samblaster} --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
	| ${samtools} view -S -b - \
	1> ${sample}-org.bam \
	2> log_path/${sample}-org.log
```


##2.2 Bam_sort
```shell
mkdir tmp_path
mkdir log_path

${sambamba} sort \
	-t ${cpus} \
	-m 10G \
	--tmpdir tmp_path/ \
	-o ${sample}-sort.bam \
	${sample}-org.bam \
	2> log_path/${sample}-sort.log


${sambamba} index \
	-t ${cpus} \
	${sample}-sort.bam \
	2> log_path/${sample}-index.log

```


##2.3 BaseRecalibrator
```java
mkdir tmp_path
mkdir log_path

java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
-XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
-Xloggc:log_path/${sample}_gc_log.BaseRecalibrator.log \
-Xms5G -Xmx15G -Djava.io.tmpdir=tmp_path \
-jar GenomeAnalysisTK.jar -T BaseRecalibrator \
-I ${sample}-sort.bam \
-o ${sample}-sort.grp \
-R GRCh37.fa \
--knownSites dbsnp-147.vcf.gz \
--knownSites Mills_and_1000G_gold_standard.indels.vcf.gz \
-nct 4 \
--read_filter BadCigar --read_filter NotPrimaryAlignment \
2> log_path/${sample}-BQSR.log
```

##2.4 PrintReads
```java
mkdir tmp_path
mkdir log_path

java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
-XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
-Xloggc:log_path/${sample}_gc_log.PrintReads.log \
-Xms5G -Xmx15G -Djava.io.tmpdir=tmp_path \
-jar GenomeAnalysisTK.jar -T PrintReads \
-I ${sample}-sort.bam \
-o ${sample}-recal-analysis.bam \
-BQSR ${sample}-sort.grp \
-EOQ \
-R GRCh37.fa \
2> log_path/${sample}-recal-analysis.log
```

