### Software versions
```
bwa 0.7.17
deeptools 3.3.0
picard 2.20.3 
trinity 2.1.1
blast 2.9.0
seqkit 0.12.0
sambamba 0.7.0
salmon 0.8.1
igblast 1.14.0
jellyfish 2.2.10
java 1.6.0
samtools 1.2
htslib 1.2.1
bedtools 2.29.0
R 4.1.1

```

### Commands

#### Alignment of short reads

```
# index genome
bwa index reference_genome.fa

# align reads, convert sam to bam file and sort by coordinates
bwa mem -t 24 bwa_genome_index file_R1.fastq file_R2.fastq | sambamba view -S -h -f bam -t 24 /dev/stdin | sambamba sort -t 24 --tmpdir=/path/ /dev/stdin -o /outpath/file.sorted.bam

# mark duplicated reads
java -jar -Xmx10G /home/pubseq/BioSw/picard/picard-tools-1.52/MarkDuplicates.jar I=file.sorted.bam O=file.sorted.dups_marked.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true

```
#### ChIP-seq track normalization

```
#deeptools
bamCoverage -b file.bam -o /outpath/file.bam.bw -of bigwig -bs 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --extendReads --ignoreDuplicates -p max/2

```