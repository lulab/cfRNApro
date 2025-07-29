# cfRNApro

## Abstract

Transcriptional dysregulation in disease is well-studied, yet how this dysregulation is reflected in cell-free RNA (cfRNA) profiles, remains largely unexplored, hindering its potential for precision medicine. Plus, A critical challenge in cfRNA analysis is the pervasive influence of technical batch effects. To address these limitations, we developed cfRNA PROcessing score (cfRNApro). This metric, derived from the statistical model of cfRNA fragmentation patterns, captures post-transcriptional regulatory signatures that show strong evidence of reflecting aberrant mRNA splicing in cancer, providing information orthogonal to conventional gene expression analyses. cfRNApro provides information orthogonal to conventional gene expression analyses and more robust to technical batch effects. The application of cfRNApro can significantly improve diagnosis and prognosis in different types of cancer.

There are three options for cfRNApro calculation: KL divergence, KS statistics and intron ratio. In our study, we prefer KL divergence because its overall performance in batch-robustness and classification, of platelet and plasma datasets, are better. However, the proformance can vary when it comes to different application situation. Thus, you may choose the option best for your research.

## Usage

