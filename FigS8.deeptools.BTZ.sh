wig_dir=/home/users/astar/gis/yaoxs/scratch/PBRM1/chip/eap/corrected_wig
deeptools_dir=/scratch/users/astar/gis/yaoxs/PBRM1/chip/deeptools

RELA_WT_DMSO=RCC591_corrected.wig.bw
RELA_WT_BTZ=RCC592_corrected.wig.bw
RELA_KO_DMSO=RCC593_corrected.wig.bw
RELA_KO_BTZ=RCC594_corrected.wig.bw

BRG1_WT_DMSO=RCC599_corrected.wig.bw
BRG1_WT_BTZ=RCC600_corrected.wig.bw
BRG1_KO_DMSO=RCC601_corrected.wig.bw
BRG1_KO_BTZ=RCC602_corrected.wig.bw

RELA_loss=/home/users/astar/gis/yaoxs/scratch/PBRM1/chip/diffbind/project213.diffbind.786O.RELA.KO.BTZ.lost.DBA_EDGER.bed

RELA_gain=/home/users/astar/gis/yaoxs/scratch/PBRM1/chip/diffbind/project213.diffbind.786O.RELA.KO.BTZ.gain.DBA_EDGER.bed

computeMatrix reference-point --referencePoint center -S $wig_dir/$RELA_WT_DMSO $wig_dir/$RELA_WT_BTZ  $wig_dir/$RELA_KO_DMSO $wig_dir/$RELA_KO_BTZ $wig_dir/$BRG1_WT_DMSO $wig_dir/$BRG1_WT_BTZ $wig_dir/$BRG1_KO_DMSO $wig_dir/$BRG1_KO_BTZ -b 1000 -a 1000 -R $RELA_loss $RELA_gain --skipZeros -o $deeptools_dir/BTZ_RELA.gz -p max/2

plotHeatmap -m $deeptools_dir/BTZ_RELA.gz -out $deeptools_dir/BRG1_RELA.pdf --colorMap Greys Blues Greys Blues --samplesLabel RELA_WT_DMSO RELA_WT_BTZ RELA_KO_DMSO RELA_KO_BTZ  SMARCA4_WT_DMSO SMARCA4_WT_BTZ SMARCA4_KO_DMSO SMARCA4_KO_BTZ --zMin 0 --yMin 0

