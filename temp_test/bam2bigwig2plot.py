BAM_COVERAGE = '/home/zhaohuanan/zhaohn_HD/miniconda3/bin/bamCoverage'
COMPUTE_MATRIX = '/home/zhaohuanan/zhaohn_HD/miniconda3/bin/computeMatrix'
PLOT_HEATMAP = '/home/zhaohuanan/zhaohn_HD/miniconda3/bin/plotHeatmap'




# step 1: DetectSeq bam to bigwig CPM and RPKM
# SAMPLE = [
#     'EMX1-All-Input_rep1',
#     'EMX1-All-Input_rep2',
#     'EMX1-All-PD_rep1',
#     'EMX1-All-PD_rep2',
#     'HN-HEK4-All-Input_rep1',
#     'HN-HEK4-All-Input_rep2',
#     'HN-HEK4-All-PD_rep1',
#     'HN-HEK4-All-PD_rep2',
#     'VEGFA-All-Input_rep1',
#     'VEGFA-All-Input_rep2',
#     'VEGFA-All-PD_rep1',
#     'VEGFA-All-PD_rep2'
# ]

# rule all:
#     input:
#         expand("../cbe_bam/293T-bat_{sample}_hg38.MAPQ20.bam", sample=SAMPLE),
#         expand("../bigwig/DetectSeq_{sample}.RPKM.bw", sample=SAMPLE),
#         expand("../bigwig/DetectSeq_{sample}.CPM.bw", sample=SAMPLE),

# rule bam2bigwig_RPKM:
#     input:
#         "../cbe_bam/293T-bat_{sample}_hg38.MAPQ20.bam"
#     output:
#         "../bigwig/DetectSeq_{sample}.RPKM.bw"
#     shell:
#         """
#         {BAM_COVERAGE} --bam {input} -o {output} -of bigwig --scaleFactor 1 --binSize 10 -p 24 --normalizeUsing RPKM
#         """
# rule bam2bigwig_CPM:
#     input:
#         "../cbe_bam/293T-bat_{sample}_hg38.MAPQ20.bam"
#     output:
#         "../bigwig/DetectSeq_{sample}.CPM.bw"
#     shell:
#         """
#         {BAM_COVERAGE} --bam {input} -o {output} -of bigwig --scaleFactor 1 --binSize 10 -p 24 --normalizeUsing CPM
#         """




# step 2: DigenomeSeq bam to bigwig CPM and RPKM



# SAMPLE = [
#     "EMX1",
#     "HEK4"
# ]

# rule all:
#     input:
#         expand("../bam/20201117-293T-WGS_.DigenomeSeq_{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam", sample=SAMPLE),
#         expand("../bigwig/DigenomeSeq_{sample}.RPKM.bw", sample=SAMPLE),
#         expand("../bigwig/DigenomeSeq_{sample}.CPM.bw", sample=SAMPLE),

# rule bam2bigwig_RPKM:
#     input:
#         "../bam/20201117-293T-WGS_.DigenomeSeq_{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam"
#     output:
#         "../bigwig/DigenomeSeq_{sample}.RPKM.bw"
#     shell:
#         """
#         {BAM_COVERAGE} --bam {input} -o {output} -of bigwig --scaleFactor 1 --binSize 10 -p 24 --normalizeUsing RPKM
#         """
# rule bam2bigwig_CPM:
#     input:
#         "../bam/20201117-293T-WGS_.DigenomeSeq_{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam"
#     output:
#         "../bigwig/DigenomeSeq_{sample}.CPM.bw"
#     shell:
#         """
#         {BAM_COVERAGE} --bam {input} -o {output} -of bigwig --scaleFactor 1 --binSize 10 -p 24 --normalizeUsing CPM
#         """







# step 3: bigwig compute matrix(CPM and RPKM)






# step 4: matrix matrix plot heatmap(CPM and RPKM)




# 293T plot 
{PLOT_HEATMAP} --perGroup \
-m deeptools_result/20200620-293T-ATACSeq_and_DNaseSeq-DetectSeq_regions.matrix.mat.gz \
-out deeptools_result/20200620-293T-ATACSeq_and_DNaseSeq.heatmap.pdf \
--regionsLabel \
DetectSeq.Out \
DetectSeq.Target \
DetectSeq.NoOutOrTarget \
DetectSeq.All \
Genome.random \
CasOFFinder.EMX1 \
CasOFFinder.VEGFA \
CasOFFinder.HEK4 &