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
# if you do not trim, just make a soft link like this

ln -s fastq fix.fastq
```
## Run
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
# see https://github.com/hermanzhaozzzz/slurm for more details
snakemake --profile slurm -pr -j 1 -s Snakefile.py -n
snakemake --profile slurm -pr -j 1 -s Snakefile.py




# use this command to see the running information
tail -f ../*/*.log
```

![](https://tva1.sinaimg.cn/large/0081Kckwly1gl0dky3r6gj30no0abq6c.jpg)

## Update logs:
- 2020-11-24: 
    - update to a normal format and remove log files
