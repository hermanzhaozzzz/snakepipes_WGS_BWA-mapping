#!/bin/bash
# properties = {"type": "single", "rule": "TruSeq_cutadapt", "local": false, "input": ["../fix.fastq/20201124-293T-DetectSeq_.CBE-V-PD7_R1.fastq.gz", "../fix.fastq/20201124-293T-DetectSeq_.CBE-V-PD7_R2.fastq.gz"], "output": ["../fix.fastq/20201124-293T-DetectSeq_.CBE-V-PD7_R1_cutadapt.fq.gz", "../fix.fastq/20201124-293T-DetectSeq_.CBE-V-PD7_R2_cutadapt.fq.gz"], "wildcards": {"date": "20201124", "cell": "293T", "seq": "DetectSeq", "sample": "CBE-V-PD7"}, "params": {}, "log": ["../fix.fastq/20201124-293T-DetectSeq_.CBE-V-PD7_cutadapt.log"], "threads": 1, "resources": {}, "jobid": 1, "cluster": {}}
 cd /gpfs/user/zhaohuanan/3.project/23.2020_11_1120_CBE-HEK4-DigenomeOnly-GuideOnly_and_LZC-CBE-DetectSeqTest-VEGFA/CBE-DetectSeqTest/snakepipes_WGS_BWA-mapping && \
PATH='/home/zhaohuanan/zhaohn_HD/miniconda3/bin':$PATH /home/zhaohuanan/zhaohn_HD/miniconda3/bin/python3.8 \
-m snakemake ../fix.fastq/20201124-293T-DetectSeq_.CBE-V-PD7_R1_cutadapt.fq.gz --snakefile /gpfs/user/zhaohuanan/3.project/23.2020_11_1120_CBE-HEK4-DigenomeOnly-GuideOnly_and_LZC-CBE-DetectSeqTest-VEGFA/CBE-DetectSeqTest/snakepipes_WGS_BWA-mapping/Snakefile.py \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/user/zhaohuanan/3.project/23.2020_11_1120_CBE-HEK4-DigenomeOnly-GuideOnly_and_LZC-CBE-DetectSeqTest-VEGFA/CBE-DetectSeqTest/snakepipes_WGS_BWA-mapping/.snakemake/tmp.ja_456sq ../fix.fastq/20201124-293T-DetectSeq_.CBE-V-PD7_R1.fastq.gz ../fix.fastq/20201124-293T-DetectSeq_.CBE-V-PD7_R2.fastq.gz --latency-wait 60 \
 --attempt 1 --force-use-threads --scheduler ilp \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules TruSeq_cutadapt --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

