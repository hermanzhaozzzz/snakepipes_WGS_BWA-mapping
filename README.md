---
Author: Herman Zhao
Email: hermanzhaozzzz@gmail.com
---


# snakepipes_WGS_BAW-mapping
## Pretreat
1. do fastqc and multiqc to check the quality of sequencing; see https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc 
```
git clone git@github.com:hermanzhaozzzz/snakepipes_fastqc-multiqc.git
```
2. if necessary , run trim protocol to make a better raw sequencing file
```
# the trim protocol will form a fix.fastq folder
# if you do not trim, just make a soft link is okay

ln -s fastq fix.fastq
```
3. fill the SAMPLE list in the snakefile
4. check the <genome, bwa index file path> in the snakefile
5. run alignment
```
cd snakepipes_WGS_BWA-mapping
# test
snakemake -pr -j 1 -s Snakefile.py -n
# run
snakemake -pr -j 1 -s Snakefile.py


# if use slurm
snakemake --profile slurm -pr -j 1 -s Snakefile.py -n
snakemake --profile slurm -pr -j 1 -s Snakefile.py




# use this command to see the running information
tail -f ../*/*.log
```

![](https://tva1.sinaimg.cn/large/0081Kckwly1gl0dd0zb0ej312f0gd41j.jpg)

## Update logs:
- 2020-11-24: 
    - update to a normal format and remove log files