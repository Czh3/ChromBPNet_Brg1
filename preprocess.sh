# generate ATAC peaks
bedtools slop -i mm10.blacklist.bed.gz -g ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.chrom.size -b 1057 > temp.bed
bedtools intersect -v -a ~/proj/Brg1_dTag/Brg1_project/ATAC/peaks/WT_ATAC/WT_ATAC_peaks.narrowPeak -b temp.bed  > WT_ATAC_peaks_no_blacklist.bed


# Split the genome into test and validation
mkdir splits
chrombpnet prep splits -c ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.chrom.size -tcr chr1 chr5 chr10 chr15 -vcr chr4 chr12 -op splits/fold_0

# Generate non-peaks (background regions)
chrombpnet prep nonpeaks -g ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -p WT_ATAC_peaks_no_blacklist.bed -c ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.chrom.size -npr 1 -fl splits/fold_0.json -br mm10.blacklist.bed.gz -o output


# Bias model training
chrombpnet bias pipeline \
  -ibam /home/zhangc/proj/Brg1_dTag/Brg1_project/ATAC/BW/RPKM/merge/WT_ATAC.bam \
  -d "ATAC" \
  -g ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
  -c ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.chrom.size \
  -p WT_ATAC_peaks_no_blacklist.bed \
  -n output_negatives.bed \
  -fl splits/fold_0.json \
  -b 0.5 \
  -o bias_model/ \
  -fp ESC_WT




########## All Brg1-binding peaks
### Ctrl Sample
# ChromBPNet training
chrombpnet pipeline \
  -ibam /home/zhangc/proj/Brg1_dTag/Brg1_project/ATAC/BW/RPKM/merge/Brg1_Ctrl.bam \
  -d "ATAC" \
  -g ~/genome/ Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
  -c ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.chrom.size \
  -p WT_ATAC_peaks_no_blacklist.bed \
  -n output_negatives.bed \
  -fl splits/fold_0.json \
  -b bias_model/models/ESC_WT_bias.h5 \
  -o Brg1_Ctrl_model

# Generate prediction bigwigs
chrombpnet contribs_bw -m Brg1_Ctrl_model/models/chrombpnet_nobias.h5 -r WT_ATAC_peaks_no_blacklist.bed        \
        -g ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa  \
        -c ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.chrom.size \
        -op Brg1_Ctrl.contribs_bw

# Denovo motif discovery
modisco motifs -i Brg1_Ctrl.contribs_bw.profile_scores.h5 -n 10000 -o  Brg1_Ctrl_motif.h5
modisco report -i  Brg1_Ctrl_motif.h5 -o  Brg1_Ctrl_report/ -s Brg1_Ctrl_report/ -m motif.meme.txt


### 100n
# ChromBPNet training
chrombpnet pipeline \
  -ibam /home/zhangc/proj/Brg1_dTag/Brg1_project/ATAC/BW/RPKM/merge/Brg1KI_24hr_dTag13_1_10_7.bam \
  -d "ATAC" \
  -g ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
  -c ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.chrom.size \
  -p WT_ATAC_peaks_no_blacklist.bed \
  -n output_negatives.bed \
  -fl splits/fold_0.json \
  -b bias_model/models/ESC_WT_bias.h5 \
  -o Brg1_100n_model

# Generate prediction bigwigs
chrombpnet contribs_bw -m Brg1_100n_model/models/chrombpnet_nobias.h5 -r WT_ATAC_peaks_no_blacklist.bed        \
        -g ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa  \
        -c ~/genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.chrom.size \
        -op Brg1_100n.contribs_bw

# Denovo motif discovery
modisco motifs -i Brg1_100n.contribs_bw.profile_scores.h5 -n 10000 -o  Brg1_100n_motif.h5
modisco report -i Brg1_100n_motif.h5 -o  Brg1_100n_report/ -s Brg1_100n_report/ -m motif.meme.txt


