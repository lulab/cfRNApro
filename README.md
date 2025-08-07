# cfRNApro

## Abstract

Transcriptional dysregulation in disease is well-studied, yet how this dysregulation is reflected in cell-free RNA (cfRNA) profiles remains largely unexplored, hindering its potential for precision medicine. Plus, A critical challenge in cfRNA analysis is the pervasive influence of technical batch effects. To address these limitations, we developed the cfRNA PROcessing score (cfRNApro). This metric, derived from the statistical model of cfRNA fragmentation patterns, captures post-transcriptional regulatory signatures that show strong evidence of reflecting aberrant mRNA splicing in cancer, providing information orthogonal to conventional gene expression analyses. cfRNApro provides information orthogonal to traditional expression of gene analyses and is more robust to technical batch effects. The application of cfRNApro can significantly improve diagnosis and prognosis in different types of cancer.

There are three options for cfRNApro calculation: KL divergence, KS statistics, and intron ratio. In our study, we prefer KL divergence because its overall performance in batch-robustness and classification of platelet and plasma datasets are better. However, the performance can vary when it comes to different application situations. Thus, you may choose the option best for your research.

## Usage

The scripts applied multiprocessing on chromosomes for acceleration. While KL and KS are independent of the reference file, IR requires a BED file describing the exon regions.

### Environment
Python 3.10.12
numpy 1.26.4
pandas

### Command

``` python
python multiprocess_KL.py --bed [bedgraph file] --out [output file]
### RETURN: 4 columns: gene names, gene length, gene abundance, and kl.

python multiprocess_KS.py --bed [bedgraph file] --out [output file]
### RETURN: 5 columns: gene names, gene length, gene abundance, max_dev (the actual statistics we used, representing the max 3' bias), and ks_statistics.

python multiprocess_IR.py --ref [reference BED file describing exon regions] --bed [bedgraph file] --out [output file]
### RETURN: 2 columns: gene names and intron_ratio.
```
