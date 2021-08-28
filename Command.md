#### alignment of short reads

```
bwa mem 

```
#### ChIP-seq track normalization

```
#deeptools
bamCoverage -b file.bam -o /outpath/file.bam.bw -of bigwig -bs 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --extendReads --ignoreDuplicates -p max/2

```