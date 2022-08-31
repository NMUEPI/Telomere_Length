# Telomere length Estimation

*Based on the whole-genome sequencing data, telseq(https://github.com/zd1/telseq) and Qmotif (https://github.com/AdamaJava/adamajava/tree/master/qmotif) were used to estimae the telomere length.*

## Telseq 

```shell
telseq -u -k 12 -r 150 -o /public/home/zj2020/telseq_res/telseq/${ID}.tel  ${ID}.bam
```


## Qmotif

```java
java -Xmx15g -jar /qmotif-1.0.jar \
     -n 4 \
     --i compare.ini  \
     --bam ${sample}-recal.bam \
     --bai ${sample}-recal.bam.bai \
     -o ${sample}.xml \
     -o ${sample}-tel.bam \
     --log ${sample}.log \
     --loglevel INFO
```

### *Configure file used for qmotif analysis*

```shell
[PARAMS]
stage1_motif_string=TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG
stage2_motif_regex=(TTAGGG|ATAGGG|CTAGGG|GTAGGG|TAAGGG|TCAGGG|TGAGGG|TTCGGG|TTGGGG|TTTGGG|TTAAGG|TTACGG|TTATGG|TTAGAG|TTAGCG|TTAGTG|TTAGGA|TTAGGC|TTAGGT|CCCTAA|ACCTAA|GCCTAA|TCCTAA|CACTAA|CGCTAA|CTCTAA|CCATAA|CCGTAA|CCTTAA|CCCAAA|CCCCAA|CCCGAA|CCCTCA|CCCTGA|CCCTTA|CCCTAC|CCCTAG|CCCTAT)
revcomp=true
window_size=10000
cutoff_size=5

[INCLUDES]
; name, regions (sequence:start-stop)
chr1p	1:10001-12464
chr1q	1:249237907-249240620
chr2p	2:10001-12592
chr2q	2:243187373-243189372
chr2xA	2:243150480-243154648
chr3p	3:60001-62000
chr3q	3:197960430-197962429
chr3xB	3:197897576-197903397
chr4p	4:10001-12193
chr4q	4:191041613-191044275
chr5p	5:10001-13806
chr5q	5:180903260-180905259
chr6p	6:60001-62000
chr6q	6:171053067-171055066
chr7p	7:10001-12238
chr7q	7:159126558-159128662
chr8p	8:10001-12000
chr8q	8:146302022-146304021
chr9p	9:10001-12359
chr9q	9:141151431-141153430
chr10p	10:60001-62000
chr10q	10:135522469-135524746
chr11p	11:60001-62000
chr11q	11:134944458-134946515
chr12p	12:60001-62000
chr12q	12:133839458-133841894
chr12xC	12:93158-97735
chr13p	13:19020001-19022000
chr13q	13:115107878-115109877
chr14p	14:19020001-19022000
chr14q	14:107287540-107289539
chr15p	15:20000001-20002000
chr15q	15:102518969-102521391
chr16p	16:60001-62033
chr16q	16:90292753-90294752
chr17p	17:1-2000
chr17q	17:81193211-81195210
chr18p	18:10001-12621
chr18q	18:78014226-78017247
chr19p	19:60001-62000
chr19q	19:59116822-59118982
chr20p	20:60001-62000
chr20q	20:62963520-62965519
chr21p	21:9411194-9413193
chr21q	21:48117788-48119894
chr22p	22:16050001-16052000
chr22q	22:51242566-51244565
chrXp	X:60001-62033
chrXq	X:155257733-155260559
chrYp	Y:10001-12033
chrYq	Y:59360739-59363565
;..

[EXCLUDES]
; regions (sequence:start-stop)

```
