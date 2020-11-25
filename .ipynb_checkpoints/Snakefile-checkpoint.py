# _*_ coding: UTF-8 _*_
import os
########################################################################
# ZHAO Huanan
# Whole Genome Sequencing Former Analysis
# [Form a bam file and do necessary filtering]
########################################################################

# --------------------------------------------------------------->>>>>>>
# pipeline
# --------------------------------------------------------------->>>>>>>
# 1. cutadapt
# 2. bwa mem mapping
# 3. bam MAPQ filter

# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
SAMPLES = [
"CBE-V-PD7"
]

# multiprocess
THREADS = '24'
# Date
DATE = '20201124'
# CELL TYPE
CELL = "293T"
# SEQ TYPE
SEQ = 'DetectSeq'









# --------------------------------------------------------------->>>>>>>
# genome, bwa index file path
# --------------------------------------------------------------->>>>>>>
# genome
HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
# bwa index 
# to form bwa index, run `bwa index hg38_only_chromosome.fa`
HG38_FA_BWA_INDEX = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"



# --------------------------------------------------------------->>>>>>>
# software
# --------------------------------------------------------------->>>>>>>
# make sure the cutadapt, bwa, samtools and sambamba are in you PATH
# get the application path
with os.popen("which cutadapt") as path:
    CUTADAPT = path.read().strip()
    print('PATH cutadapt:', CUTADAPT)
with os.popen("which bwa") as path:
    BWA = path.read().strip()
    print('PATH bwa:', BWA)
with os.popen("which samtools") as path:
    SAMTOOLS = path.read().strip()
    print('PATH samtools:', SAMTOOLS)
with os.popen("which sambamba") as path:
    SAMBAMBA = path.read().strip()
    print('PATH sambamba:', SAMBAMBA)
    
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
        {HG38_FA_BWA_INDEX} \
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
    log:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.bam.log"
    shell:
        """
        {SAMTOOLS} view -@ {THREADS} -Sbh \
        {input} > {output} > {log} 2>&1
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
        {SAMTOOLS} addreplacerg -r {params.tag} -@ {THREADS} -O BAM -o {output} \
        --reference {HG38_FA} {input}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools sort by position
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_sort_by_position:
    input:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned.out.fix_RG.bam"
    output:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.bam"
    log:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.bam.log"
    shell:
        """
        {SAMTOOLS} sort \
        -O BAM \
        -o {output} \
        -T {output}.temp \
        -@ {THREADS}  \
        {input} \
        > {log} 2>&1
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
# sambamba rmdup and build bam index
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_index:
    input:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20.bam"
    output:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam",
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam.bai"
    log:
        "../bam/{date}-{cell}-{seq}_.{sample}_Aligned_sortp.mapped_F1804_MAPQ20_rmdup.bam.log",
    shell:
        """
        {SAMBAMBA} markdup --remove-duplicates --nthreads= {THREADS} --show-progress  --sort-buffer-size 8192 \
        {input} {output[0]} > {log} 2>&1
        """
        