#A promoter region is defined as H3K27ac positive, H3K4me3 positive and within 2.5kb +/- TSS where TSS - GENCODE transcripts
#A enhancer region is defined as H3K27ac positive, H3K4me1 positive and non-overlapping with promoter 

out=/mnt/projects/yaoxs/cancer_chipseqxs/PBRM1/intersect
TSS_2000_GENCODE=/mnt/projects/yaoxs/cancer_chipseqxs/kidney/others/UCSC.hg19.knownGene.tss.4kb.bed
peak_dir=/mnt/projects/tanbop/intergrated_genomics/rcc/chip/MACS2/

# histone peak regions
H3K27ac_WT=project197.mem/RCC431.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC407.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed
H3K27ac_KO=project197.mem/RCC432.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC429.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed
H3K4me3_WT=project197.mem/RCC435.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC407.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed
H3K4me3_KO=project197.mem/RCC436.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC429.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed
H3K4me1_WT=project208.mem/RCC542.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC543.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed
H3K4me1_KO=project208.mem/RCC547.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC548.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed


# SWI/SNF peak regions

#PBRM1
PBRM1_1=project179.mem/RCC346.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC341.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed
PBRM1_2=project194.mem/RCC406.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC407.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed
PBRM1_3=project194.mem/RCC408.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC409.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed


#ARID1A
ARID1A_1=project211.mem/RCC575.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC543.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed
ARID1A_2=project211.mem/RCC576.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC548.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed


#ARID2
ARID2_1=project188.mem/RCC366_peaks.bed
ARID2_2=project211.mem/RCC578.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC543.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed

#BRD7
BRD7_1=project209.mem/RCC562.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC543.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed
BRD7_2=project209.mem/RCC570.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC548.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed

#SMARCA4
SMARCA4_1=project188.mem/RCC364_peaks.bed
SMARCA4_2=project211.mem/RCC577.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC543.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed


#merge H3K27ac
#merge H3K4me3
#merge H3K4me1


K27ac=K27ac.merge.all.bed
K4me3=K4me3.merge.all.bed
K4me1=K4me1.merge.all.bed

cat $peak_dir/$H3K27ac_WT $peak_dir/$H3K27ac_KO |sortBed|mergeBed > $out/$K27ac
cat $peak_dir/$H3K4me3_WT $peak_dir/$H3K4me3_KO |sortBed|mergeBed > $out/$K4me3
cat $peak_dir/$H3K4me1_WT $peak_dir/$H3K4me1_KO |sortBed|mergeBed > $out/$K4me1



#promoters
bedtools intersect -a $out/"$K27ac" -b $out/"$K4me3" > $out/K27ac.K4me3.bed
bedtools intersect -u -a $out/K27ac.K4me3.bed -b $TSS_2000_GENCODE |uniq > $out/K27ac.K4me3.TSS_2000_GENCODE.bed
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "prom-"(NR+1000000)}' $out/K27ac.K4me3.TSS_2000_GENCODE.bed > $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed


#enhancers
bedtools intersect -a $out/"$K27ac" -b $out/"$K4me1" > $out/K27ac.K4me1.bed
bedtools intersect -v -a $out/K27ac.K4me1.bed -b $out/K27ac.K4me3.TSS_2000_GENCODE.bed |uniq > $out/K27ac.K4me1.noprom.bed
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "enh-"(NR+1000000)}' $out/K27ac.K4me1.noprom.bed > $out/K27ac.K4me1.noprom.ID.bed

 
#perform overlap with BRG1 peaks with promoters and enhancers



bedtools intersect -wa -a $peak_dir/$PBRM1_1 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/PBRM1_1.promoter.bed
bedtools intersect -wa -a $peak_dir/$PBRM1_1 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >	 $out/PBRM1_1.enhancer.bed

bedtools intersect -wa -a $peak_dir/$PBRM1_2 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/PBRM1_2.promoter.bed
bedtools intersect -wa -a $peak_dir/$PBRM1_2 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >     $out/PBRM1_2.enhancer.bed

bedtools intersect -wa -a $peak_dir/$PBRM1_3 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/PBRM1_3.promoter.bed
bedtools intersect -wa -a $peak_dir/$PBRM1_3 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >     $out/PBRM1_3.enhancer.bed




bedtools intersect -wa -a $peak_dir/$ARID1A_1 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/ARID1A_1.promoter.bed
bedtools intersect -wa -a $peak_dir/$ARID1A_1 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >     $out/ARID1A_1.enhancer.bed

bedtools intersect -wa -a $peak_dir/$ARID1A_2 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/ARID1A_2.promoter.bed
bedtools intersect -wa -a $peak_dir/$ARID1A_2 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >     $out/ARID1A_2.enhancer.bed


bedtools intersect -wa -a $peak_dir/$ARID2_1 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/ARID2_1.promoter.bed
bedtools intersect -wa -a $peak_dir/$ARID2_1 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >     $out/ARID2_1.enhancer.bed

bedtools intersect -wa -a $peak_dir/$ARID2_2 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/ARID2_2.promoter.bed
bedtools intersect -wa -a $peak_dir/$ARID2_2 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >     $out/ARID2_2.enhancer.bed

bedtools intersect -wa -a $peak_dir/$BRD7_1 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/BRD7_1.promoter.bed
bedtools intersect -wa -a $peak_dir/$BRD7_1 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >     $out/BRD7_1.enhancer.bed

bedtools intersect -wa -a $peak_dir/$BRD7_2 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/BRD7_2.promoter.bed
bedtools intersect -wa -a $peak_dir/$BRD7_2 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >     $out/BRD7_2.enhancer.bed

bedtools intersect -wa -a $peak_dir/$SMARCA4_1 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/SMARCA4_1.promoter.bed
bedtools intersect -wa -a $peak_dir/$SMARCA4_1 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >     $out/SMARCA4_1.enhancer.bed

bedtools intersect -wa -a $peak_dir/$SMARCA4_2 -b $out/K27ac.K4me3.TSS_2000_GENCODE.ID.bed |uniq >  $out/SMARCA4_2.promoter.bed
bedtools intersect -wa -a $peak_dir/$SMARCA4_2 -b $out/K27ac.K4me1.noprom.ID.bed |uniq >     $out/SMARCA4_2.enhancer.bed

