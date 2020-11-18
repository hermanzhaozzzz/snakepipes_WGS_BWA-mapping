#!/bin/bash
# properties = {"type": "single", "rule": "pick_mapped_reads", "local": false, "input": ["../bam/20201117-293T-WGS_.DigenomeSeq_HEK4_Aligned_sortp.bam"], "output": ["../bam/20201117-293T-WGS_.DigenomeSeq_HEK4_Aligned_sortp.mapped_F1804_MAPQ20.bam"], "wildcards": {"date": "20201117", "cell": "293T", "seq": "WGS", "sample": "DigenomeSeq_HEK4"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 20, "cluster": {}}
 cd /gpfs/user/zhaohuanan/3.project/22.2020-11_1112_JSK_DigenomeSeq_WGS/snakepipes_WGS-BWAmapping && \
PATH='/home/zhaohuanan/zhaohn_HD/miniconda3/bin':$PATH /home/zhaohuanan/zhaohn_HD/miniconda3/bin/python3.8 \
-m snakemake ../bam/20201117-293T-WGS_.DigenomeSeq_HEK4_Aligned_sortp.mapped_F1804_MAPQ20.bam --snakefile /gpfs/user/zhaohuanan/3.project/22.2020-11_1112_JSK_DigenomeSeq_WGS/snakepipes_WGS-BWAmapping/Snakefile.py \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/user/zhaohuanan/3.project/22.2020-11_1112_JSK_DigenomeSeq_WGS/snakepipes_WGS-BWAmapping/.snakemake/tmp.q5mwqwtn ../bam/20201117-293T-WGS_.DigenomeSeq_HEK4_Aligned_sortp.bam --latency-wait 60 \
 --attempt 1 --force-use-threads --scheduler ilp \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules pick_mapped_reads --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

