#!/bin/bash
# properties = {"type": "single", "rule": "pick_unmapped_reads1", "local": false, "input": ["../bam/20201117-293T-WGS_.DigenomeSeq_EMX1_Aligned_sortp.bam"], "output": ["../bam/20201117-293T-WGS_.DigenomeSeq_EMX1_Aligned_sortp.unmapped.bam"], "wildcards": {"date": "20201117", "cell": "293T", "seq": "WGS", "sample": "DigenomeSeq_EMX1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 13, "cluster": {}}
 cd /gpfs/user/zhaohuanan/3.project/22.2020-11_1112_JSK_DigenomeSeq_WGS/snakepipes_WGS-BWAmapping && \
PATH='/home/zhaohuanan/zhaohn_HD/miniconda3/bin':$PATH /home/zhaohuanan/zhaohn_HD/miniconda3/bin/python3.8 \
-m snakemake ../bam/20201117-293T-WGS_.DigenomeSeq_EMX1_Aligned_sortp.unmapped.bam --snakefile /gpfs/user/zhaohuanan/3.project/22.2020-11_1112_JSK_DigenomeSeq_WGS/snakepipes_WGS-BWAmapping/Snakefile.py \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/user/zhaohuanan/3.project/22.2020-11_1112_JSK_DigenomeSeq_WGS/snakepipes_WGS-BWAmapping/.snakemake/tmp.q5mwqwtn ../bam/20201117-293T-WGS_.DigenomeSeq_EMX1_Aligned_sortp.bam --latency-wait 60 \
 --attempt 1 --force-use-threads --scheduler ilp \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules pick_unmapped_reads1 --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

