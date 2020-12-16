# VSS: Variance-stabilized units for sequencing-based genomic signals

## Variance-stabilized signals (VSS) is a signal transformation approach used for eliminating the dependency of data variance from its mean. We generate VSS for sequencing-based genomic signals by learning the empirical relationship between the mean and variance of a given signal data set and producing transformed signals that normalize for this dependence.
## Two main steps in VSS pipeline:
### 1: Training model: Identifying the mean-variance relationship
There are two different options for this step. It either uses the user provided replicates to identify the mean-variance relationship or uses the default trained model. In the latter case, user just needs to provide the untransformed signals.
### 2: Transforming signals: Calculating variance-stabilized signals
Having learned the mean-variance relationship, VSS can be generated using the variance-stabilizing transformation. 



<img src="https://github.com/faezeh-bayat/Variance-stabilized-units-for-sequencing-based-genomic-signals/blob/master/bin/VSS_general_schematic/VSS_schematic.png" width="800"/>

#####################################################################################

## Prerequisites
```
R (https://www.r-project.org/)
R packages needed
install.packages('bigmemory', 'data.table', 'argparse', 'pracma')

conda install -c macs2
conda install -c bioconda ucsc-bedclip
conda install -c bioconda ucsc-bedgraphtobigwig
conda install -c bioconda ucsc-bigwigtobedgraph

bedtools
######(https://bedtools.readthedocs.io/en/latest/content/installation.html)
```

## Installing VSS pipeline
```
git clone https://github.com/faezeh-bayat/Variance-stabilized-units-for-sequencing-based-genomic-signals.git
```

## How to run VSS pipeline
#### 1. Train the model
```
Rscript VSS.R train rep1 <bed, bedGraph, bam, bigWig> rep2.bed <bed, bedGraph, bam, bigWig> triandir
               
```

```
Rscript VSS.R train_tag \
              tag_alignment_rep1 <bam> \
              --fraglen1 <Fragment length for replicate 1>
              tag_alignment_rep2 <bam> \
              --fraglen2 <Fragment length for replicate 2> \
              --chrsz <2-col chromosome sizes file> \
              --gensz <hs, mm>
              --signal <fc, pval, both> 
              --outdir triandir
             
```
#### 1. Transform the signals
```
Rscript VSS.R transform rep1 <bed, bedGraph, bam, bigWig> outputtriandir traindir tranformdir
               
```

#### 3. Variance-stabilized signals will be saved in the tranformdir folder as "Variance_stabilized_signals.bed".


##################################################################################################

#### You can accesee the "Variance-stabilized units for sequencing-based genomic signals" manuscript in:
https://www.biorxiv.org/content/10.1101/2020.01.31.929174v2

