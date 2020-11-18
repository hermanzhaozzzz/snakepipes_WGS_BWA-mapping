# snakepipes_cutadapt-HISAT2mapping-FPKM-sortBAM


```
bedtools intersect -a DetectSeq.HEK4.bed -b DigenomeSeq.HEK4.bed -loj > DetectSeq_vs_DigenomeSeq_base-on_DetectSeq_HEK4.csv

bedtools intersect -a DigenomeSeq.HEK4.bed -b DetectSeq.HEK4.bed -loj > DetectSeq_vs_DigenomeSeq_base-on_DigenomeSeq_HEK4.csv
```

然后用setect_region.ipynb筛选出only detectseq only digenomeseq和二者share的bed文件

```
# check独立性，输出的追加列都是-1，说明无重叠
bedtools intersect -a DetectSeq_vs_DigenomeSeq_Final-DetectSeqOnly_HEK4.bed -b DetectSeq_vs_DigenomeSeq_Final-DigenomeSeqOnly_HEK4.bed -loj
```