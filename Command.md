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
macs2 2.1.1.20160309
finder 2
```

### Commands

#### Alignment of short reads

```
# index genome
bwa index reference_genome.fa

# align reads, convert sam to bam file and sort by coordinates
bwa mem -t 24 bwa_genome_index file_R1.fastq file_R2.fastq | sambamba view -S -h -f bam -t 24 /dev/stdin | sambamba sort -t 24 --tmpdir=/path/ /dev/stdin -o /outpath/file.sorted.bam

# mark duplicated reads
java -jar -Xmx10G MarkDuplicates.jar I=file.sorted.bam O=file.sorted.dups_marked.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true

# making tracks for UCSC genome browser
                                
mark=H3K36me3
for i in 1;
do
echo "track $mark"
echo "compositeTrack on"
echo "shortLabel $mark"
echo "longLabel $mark"
echo "type bigWig"
echo ""

for f in *NBC*bw *GCBC*bw *MBC*bw *PBC*bw *uCLL*bw *mCLL*bw *_CLL_*;  
do 
echo $f | awk '{gsub(".bw", "") ; print $0}' | tr -s '_' '\t' | awk '{print "        " " track", $1"_"$4"_"$5"_"$6"_"$7$8}'
echo $f | awk '{gsub(".bw", "") ; print $0}' | tr -s '_' '\t' | awk '{print "        " " shortLabel", $5"_"$6"_"$7$8}'
echo $f | awk '{gsub(".bw", "") ; print $0}' | tr -s '_' '\t' | awk '{print "        " " longLabel", $1"_"$4"_"$5"_"$6"_"$7$8}'
echo "        " "parent $mark on"
echo "        " "type bigWig"
echo "        " "visibility full"
echo "        " "maxHeightPixels 70:70:32"
echo "        " "configurable on"
echo "        " "autoScale on"
echo "        " "alwaysZero on"
echo "        " "priority 0.1"
echo "        " "bigDataUrl $f"
echo "        " "color 153,0,153"
echo "        " ""
done 
done >trackDb.txt

#link hubs
http://www.epigenomes.ca/data/CLL_rislam/H3K36me3/hub.txt 

## analysis of differential ChIP-seq regions. Differential ChIP-seq regions were called using an in-house pipeline.

# filter region: merge de regions within 100bp, min 300bp, pvalue 0.03 

cd $outdir
while read line; do 
	echo $line;
	cat $indir/$line*0.03 | sed 's/,/\t/g' | grep -v "chromosome" | awk '{print $1 "\t" $2 "\t" $2+50 "\t" $5 "\t" $8 "\t" $12 "\t" $13} ' | sort -k1,1 -k2,2n | bedtools merge -d 100 -c 4,5,6 -o mean,mean,collapse -i stdin | awk '$3-$2 >=300{print $0}' >$line.100bpMerge_300bpMin_0.03Pvalue;  
done <$file

# get cell types enriched over other
for f in *_0.03Pvalue; do 
	echo $f; 
	for cell in CLL NBC GCBC PBC MBC cll_unmutated cll_mutated del_CT wt; do 
		echo $cell; 
		echo ${f%.100bpMerge_300bpMin_0.03Pvalue}"_where_"$cell"_enriched";
		less $f  | grep $cell >${f%.100bpMerge_300bpMin_0.03Pvalue}"_where_"$cell"_enriched";
	done
done

# remove unnecessary files
ls -l *_enriched | awk '$5==0{print "rm" "\t" $9}' >remove_emptyfiles.sh
bash remove_emptyfiles.sh 
rm *_0.03Pvalue remove_emptyfiles.sh

# CLL vs 4 others
# get marks enriched (UP) in 3/4 comparisons
mkdir -p cll_bcell
for mark in H3K27ac H3K27me3 H3K4me3 H3K4me1 H3K9me3 H3K36me3; do 
	echo $mark;
	for cell in CLL NBC GCBC PBC MBC; do 
		echo $cell; 
		ls -l $mark*"where_"$cell"_enriched"; 
		bedtools multiinter -i $mark*"where_"$cell"_enriched" | awk '$4>=3{print}' | sort -k1,1 -k2,2n | bedtools merge -i stdin >./cll_bcell/$mark"_"$cell"_enriched_over_others.bed";
	done
done

# make heatmap matrix all marks
# loop: CLL + 4 bcells
cd $outdir/cll_bcell
mkdir -p matrix

for mark in H3K27ac H3K27me3 H3K4me3  H3K4me1  H3K9me3  H3K36me3; do
	for cell in CLL NBC GCBC PBC MBC; do
		ls -l $mark"_"$cell"_enriched_over_others.bed"; 
		less $mark"_"$cell"_enriched_over_others.bed" | awk '{print $1"_"$2"_"$3}' >A.intersect.value;
		for f in $bed/*$mark*bed;  do 
			outfile=$(basename $f);
			intersectBed -a <( less $mark"_"$cell"_enriched_over_others.bed" | awk '{print $1 "\t" $2 "\t" $3}') -b <(less $f | awk '{print $1 "\t" $2 "\t" $3}' ) -wao | awk '{print $1"_"$2"_"$3 "\t" $0}' | awk '!seen[$1]++' | awk '{print $8}' >$outfile.intersect.value;
		done
		((echo *intersect.value |tr ' ' '\t') && (paste *.intersect.value)) >./matrix/$mark.$cell".enriched_over4cell_0.75percent_matrix.tsv";
		rm *intersect.value
	done
done 


# remove unnecessary files
rm *intersect.value

# get marks NOT enriched (DN) in CLL in 3/4 comparisons
cd $outdir

#CLL
for mark in H3K27ac H3K27me3 H3K4me3  H3K4me1  H3K9me3  H3K36me3; do
	ls -l $mark"_CLL_"*"BC_enriched"; 
	bedtools multiinter -i $mark"_CLL_"*"BC_enriched" | awk '$4>=3{print}' | sort -k1,1 -k2,2n | bedtools merge -i stdin >./cll_bcell/$mark"_CLL_not_enriched_over_others.bed";
done

#NBC
for mark in H3K27ac H3K27me3 H3K4me3  H3K4me1  H3K9me3  H3K36me3; do
	echo $mark; 
	bedtools multiinter -i $mark"_CLL_NBC_where_CLL_enriched" $mark"_GCBC_NBC_where_GCBC_enriched" $mark"_MBC_NBC_where_MBC_enriched" $mark"_NBC_PBC_where_PBC_enriched" | awk '$4>=3{print}' | sort -k1,1 -k2,2n | bedtools merge -i stdin >./cll_bcell/$mark"_NBC_not_enriched_over_others.bed";
done

#GCBC
ll H3K9me3*GCBC*
for mark in H3K27ac H3K27me3 H3K4me3  H3K4me1  H3K9me3  H3K36me3; do
	echo $mark; 
	bedtools multiinter -i $mark"_CLL_GCBC_where_CLL_enriched" $mark"_GCBC_MBC_where_MBC_enriched" $mark"_GCBC_NBC_where_NBC_enriched" $mark"_GCBC_PBC_where_PBC_enriched" | awk '$4>=3{print}' | sort -k1,1 -k2,2n | bedtools merge -i stdin >./cll_bcell/$mark"_GCBC_not_enriched_over_others.bed";
done

#MBC
ll H3K9me3*MBC*
for mark in H3K27ac H3K27me3 H3K4me3  H3K4me1  H3K9me3  H3K36me3; do
	echo $mark; 
	bedtools multiinter -i $mark"_CLL_MBC_where_CLL_enriched" $mark"_GCBC_MBC_where_GCBC_enriched" $mark"_MBC_NBC_where_NBC_enriched" $mark"_MBC_PBC_where_PBC_enriched" | awk '$4>=3{print}' | sort -k1,1 -k2,2n | bedtools merge -i stdin >./cll_bcell/$mark"_MBC_not_enriched_over_others.bed";
done

#PBC
ll H3K9me3*PBC*
for mark in H3K27ac H3K27me3 H3K4me3  H3K4me1  H3K9me3  H3K36me3; do
	echo $mark; 
	bedtools multiinter -i $mark"_CLL_PBC_where_CLL_enriched" $mark"_GCBC_PBC_where_GCBC_enriched" $mark"_MBC_PBC_where_MBC_enriched" $mark"_NBC_PBC_where_NBC_enriched" | awk '$4>=3{print}' | sort -k1,1 -k2,2n | bedtools merge -i stdin >./cll_bcell/$mark"_PBC_not_enriched_over_others.bed";
done

# make heatmap matrix for all marks
# loop: CLL + 4 bcells
cd $outdir/cll_bcell
bed=/pathToBedfiles/

for mark in H3K27ac H3K27me3 H3K4me3  H3K4me1  H3K9me3  H3K36me3; do
		ls -l $mark"_CLL_not_enriched_over_others.bed"; 
		less $mark"_CLL_not_enriched_over_others.bed" | awk '{print $1"_"$2"_"$3}' >A.intersect.value;
		for f in $bed/*$mark*bed;  do 
			outfile=$(basename $f);
			intersectBed -a <( less $mark"_CLL_not_enriched_over_others.bed" | awk '{print $1 "\t" $2 "\t" $3}') -b <(less $f | awk '{print $1 "\t" $2 "\t" $3}' ) -wao | awk '{print $1"_"$2"_"$3 "\t" $0}' | awk '!seen[$1]++' | awk '{print $8}' >$outfile.intersect.value;
		done 
		((echo *intersect.value |tr ' ' '\t') && (paste *.intersect.value)) >./matrix/$mark"CLL_not_enriched_over4cell_0.75percent_matrix.tsv";
		rm *intersect.value
done 
```
#### ChIP-seq data analysis

```
# peak calling with findER
java -jar -Xmx25G /path/finder2.jar inputBam:$input.bam signalBam:$signal.bam outDir:$out acgtDir:/path/hg38/ACGT

# peak calling with MACS2
macs2 callpeak -B -t file.bam  -f BAM -n $line.noInput -g hs --call-summits  --outdir /path/MACS2/

# RPKM normalization of bigwig files using deeptools
bamCoverage -b file.bam -o /outpath/file.bam.bw -of bigwig -bs 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --extendReads --ignoreDuplicates -p max/2

# generate RPKM matrix using deeptools
computeMatrix reference-point --referencePoint center \
                              -S KOPTK1_RUNX1-On_H3K27ac.bam.bw \
                                 KOPTK1_RUNX1-Off_H3K27ac.bam.bw \
                              -R RUNX1_peaks.bed \
                              -a 2000 \
                              -b 2000 \
                              -p 16 \
                              --skipZeros -o matrix.gz 

# plotheatmap using deeptools
plotHeatmap -m matrix.gz \
            -out RUNX1_H3K27ac.png \
            --dpi 300 \
            --colorList '#ffeda0,blue' \
            -y 'Enrichment' \
            --heatmapWidth 5 \
            --zMin 0  --zMax 15 

```

#### RNA-seq data analysis

```
# repositioning 

# RNAseq master (in-house) 

```

#### WGBS data analysis

```
# alignment using novoalign 

# combine CpGs

# call fraction of methylation
```