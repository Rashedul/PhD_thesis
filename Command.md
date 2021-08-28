### Software versions
```

```

### Commands
#### Alignment of short reads

```
# index genome
bwa index reference_genome.fa

# align reads
bwa mem 

# sort reads 
sambamba sort -t 24 file.bam

# mark duplicated reads
java -jar -Xmx10G /home/pubseq/BioSw/picard/picard-tools-1.52/MarkDuplicates.jar I=file.sorted.bam O=file.sorted.dups_marked.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true

```
#### ChIP-seq track normalization

```
#deeptools
bamCoverage -b file.bam -o /outpath/file.bam.bw -of bigwig -bs 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --extendReads --ignoreDuplicates -p max/2

```