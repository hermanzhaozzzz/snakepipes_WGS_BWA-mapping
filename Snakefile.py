# _*_ coding: UTF-8 _*_
########################################################################
# ZHAO Huanan
# 2020-11-17
# Whole Genome Sequencing Former Analysis[Form a bam file and do necessary filtering]
######################################################################## 
# run on abyss

# before this, make sure you have done fastqc + multiqc 
# https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc
# and add trim rule to make a better trim
# --------------------------------------------------------------->>>>>>>
# pipeline
# --------------------------------------------------------------->>>>>>>
# 1. cutadapt
# 2. bwa mem mapping
# 3. bam MAPQ filter
# --------------------------------------------------------------->>>>>>>
# software
# --------------------------------------------------------------->>>>>>>
CUTADAPT = "/home/zhaohuanan/zhaohn_HD/miniconda3/bin/cutadapt"
BWA = "/home/zhaohuanan/zhaohn_HD/miniconda3/bin/bwa"
SAMTOOLS = "/home/zhaohuanan/zhaohn_HD/miniconda3/bin/samtools"
JAVA = "/home/zhaohuanan/zhaohn_HD/miniconda3/bin/java"
SAMBAMBA = "/home/zhaohuanan/zhaohn_HD/miniconda3/bin/sambamba"
PICARD = "/home/zhaohuanan/zhaohn_HD/1.apps/picard/picard.jar"
# --------------------------------------------------------------->>>>>>>
# index and files
# ------------------------------------------------------ --------->>>>>>>
# genome
HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
# bwa index
HG38_FA_DICT = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
SAMPLES = [
"DigenomeSeq_EMX1",
"DigenomeSeq_HEK4"
]
# multiprocess
THREADS = '24'
# Date
DATE = '20201117'
# CELL TYPE
CELL = "293T"
# SEQ TYPE
SEQ = 'WGS'
# --------------------------------------------------------------->>>>>>>







# --------------------------------------------------------------->>>>>>>
# start
# --------------------------------------------------------------->>>>>>>
# you better use slurm cmd to run this snk file job
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule all:
    input:
        expand("../fix.fastq/{date}-{cell}-{seq}_.{sample}_R1.fastq.gz",date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
        expand("../fix.fastq/{date}-{cell}-{seq}_.{sample}_R2.fastq.gz",date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
        expand("../fix.fastq/{date}-{cell}-{seq}_.{sample}_R1_cutadapt.fq.gz",date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
        expand("../fix.fastq/{date}-{cell}-{seq}_.{sample}_R2_cutadapt.fq.gz",date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
        expand("../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.bam",date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
        expand("../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.fix_RG.bam",date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
        expand("../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.bam",date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
        expand("../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.unmapped.fasta.gz",
               date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
        expand("../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam",
               date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
#         expand("../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.matrix",
#                date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
        expand("../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam.bai",
               date=DATE,cell=CELL,seq=SEQ,sample=SAMPLES),
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# cutadapter
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule TruSeq_cutadapt:
    input:
        "../fix.fastq/{date}-{cell}-{seq}_.{sample}_R1.fastq.gz","../fix.fastq/{date}-{cell}-{seq}_.{sample}_R2.fastq.gz"
    output:
        "../fix.fastq/{date}-{cell}-{seq}_.{sample}_R1_cutadapt.fq.gz",
        "../fix.fastq/{date}-{cell}-{seq}_.{sample}_R2_cutadapt.fq.gz"
    log:
        "../fix.fastq/{date}-{cell}-{seq}_.{sample}_cutadapt.log"
    shell:# using illumina universal adaptor
        """
        {CUTADAPT} -j {THREADS} --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 55 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# Hisat2 mapping
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BWA_mapping:
    input:
        fq1 = "../fix.fastq/{date}-{cell}-{seq}_.{sample}_R1_cutadapt.fq.gz",
        fq2 = "../fix.fastq/{date}-{cell}-{seq}_.{sample}_R2_cutadapt.fq.gz"
    output:
        temp("../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.sam")
    log:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.log"
    
    shell:
        """
        {BWA} mem \
        -t {THREADS} \
        {HG38_FA_DICT} \
        {input.fq1} \
        {input.fq2} \
        -o {output}  > {log} 2>&1
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# convert sam to bam
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule sam2bam:
    input:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.sam"
    output:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.bam"
    shell:
        """
        {SAMTOOLS} view -@ {THREADS} -Sbh \
        {input} > {output}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# add @RG tag (mostly for GATK SNP/SNV calling)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule add_RG_tag:
    input:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.bam"
    output:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.fix_RG.bam"
    params:
        tag = "'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'"
    shell:
        """
        {SAMTOOLS} addreplacerg -r {params.tag} -@ {THREADS} -O BAM -o {output} --reference {HG38_FA_DICT} {input}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools sort by position(not sort by name)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_sort_by_position:
    input:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.fix_RG.bam"
    output:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.bam"
    shell:
        """
        {SAMTOOLS} sort \
        -O BAM \
        -o {output} \
        -T {output}.temp \
        -@ {THREADS}  \
        {input}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools view -f 4 step1
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule pick_unmapped_reads1:
    input:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.bam"
    output:
        temp("../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.unmapped.bam")
    shell:
        """
        {SAMTOOLS} view -@ {THREADS} -b -f 4 \
        {input} > {output}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools view -f 4 step2
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule pick_unmapped_reads2:
    input:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.unmapped.bam"
    output:
        temp("../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.unmapped.fasta")
    params:
        awk_param = """'{print ">"$1"\\n"$10}'"""
    shell:
        """
        {SAMTOOLS} view -@ {THREADS} {input} | awk {params.awk_param} > {output}
        """     
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools view -f 4 step3
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule pick_unmapped_reads3:
    input:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.unmapped.fasta"
    output:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.unmapped.fasta.gz"
    shell:
        """
        pigz -p {THREADS} {input}
        """         
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools view -F 1804 -q 20
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule pick_mapped_reads:
    input:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.bam"
    output:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20.bam"
    shell:
        """
        {SAMTOOLS} view -@ {THREADS} -b -F 1804 -q 20  \
        {input} > {output}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# picard remove duplications
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule Picard_rmdup:
#     input:
#         "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20.bam"
#     output:
#         "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam",
#         "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.matrix"
#     log:
#         "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.log"
#     shell:
#         """
#         {JAVA} -Xms90g -Xmx90g -XX:ParallelGCThreads={THREADS} \
#         -jar {PICARD} MarkDuplicates \
#         I={input} \
#         O={output[0]} \
#         M={output[1]} \
#         ASO=coordinate \
#         REMOVE_DUPLICATES=true 2>{log}
#         """
# # ------------------------------------------------------------------------------------------>>>>>>>>>>
# # samtools build bam index
# # ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule BAM_index:
#     input:
#         "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam"
#     output:
#         "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam.bai"
#     shell:
#         """
#         {SAMTOOLS} index -@ {THREADS} \
#         {input} \
#         {output}
#         """

# ------------------------------------------------------------------------------------------>>>>>>>>>>
# sambamba rmdup and build bam index
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_index:
    input:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20.bam"
    output:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam",
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam.bai"
    shell:
        """
        {SAMBAMBA} markdup --remove-duplicates --nthreads= {THREADS} --show-progress  --sort-buffer-size 8192 \
        {input} {output[0]}
        """
        