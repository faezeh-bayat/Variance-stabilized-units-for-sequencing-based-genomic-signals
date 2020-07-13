# Variance-stabilized units for sequencing-based genomic signals

## Variance-stabilized signals (VSS) is a signal transformation approach used for eliminating the dependency of data variance from its mean. We generate VSS for sequencing-based genomic signals by learning the empirical relationship between the mean and variance of a given signal data set and producing transformed signals that normalize for this dependence.
## Two main steps in VSS pipeline:
### 1: Identifying the mean-variance relationship
There are two different options for this step. It either uses the user provided replicates to identify the mean-variance relationship or uses the default trained model. In the latter case, user just needs to provide the untransformed signals.
### 2: Calculating variance-stabilized signals
Having learned the mean-variance relationship, VSS can be generated using the variance-stabilizing transformation. 



<img src="https://github.com/faezeh-bayat/Variance-stabilized-units-for-sequencing-based-genomic-signals/blob/master/bin/VSS_general_schematic/VSS.png" width="800"/>

#####################################################################################

## Prerequisites
```
R (https://www.r-project.org/)
R packages needed
install.packages('bigmemory')
install.packages('pracma')

bedtools
######(https://bedtools.readthedocs.io/en/latest/content/installation.html)
```

## Installing VSS pipeline
```
git clone https://github.com/guanjue/S3V2_IDEAS_ESMP.git
```



## Inputs for VSS pipeline
#### 1. Users want to build their own mean-variance relationship model:
##### Two untransformed replicates should be provided for building the model. After the model is build, any untransformed signals can be provided to get the varinace-stabilized signals. All input signals can be in bedGraph, bigWig or bam format.

```
variance_stabilization_model="user_specified"
replicate1_signals_for_training_the_model="rep1.bedGraph"
replicate2_signals_for_training_the_model="rep2.bedGraph"
chromosomes_to_build_the_model="chr21"
signals_to_be_variance_stabilized="rep1.bedGraph"
chromosomes_to_be_stabilized="chr21"
```
#### 2. Users want to use the default mean-variance relationship:
##### Only one untransformed signal that needs to be variance-stabilized should be provided. This input can also be in bedGraph, bigWig or bam format.
```
variance_stabilization_model="default"
replicate1_signals_for_training_the_model="False"
replicate2_signals_for_training_the_model="False"
chromosomes_to_build_the_model="False"
signals_to_be_variance_stabilized="rep1.bedGraph"
chromosomes_to_be_stabilized="chr21"
```
## Input file format
```
chromA  chromStartA  chromEndA  dataValueA
chromB  chromStartB  chromEndB  dataValueB
```



##### It only needs one input metadata file which tells the pipeline where are input 
files. An example of the metadata is in the "metadata.for_master_peak_calls.txt" file with 4 columns (columns are separated by tab):
##### 1st column: cell type name (!!!The cell type name should not have "." in it!!!)
##### 2nd column: epigenetic feature
##### 3rd column: cell type id
##### 4th column: absolute path to the IP bigwig files
##### 5th column: absolute path to the CONTROL bigwig files (If there is no control signal track, this column can be leave as empty)
```
>>> head metadata.for_master_peak_calls.txt
Cell1	ATAC	01	/storage/home/gzx103/scratch/S3V2norm_compare/hg38_cCREs/bw/c10_chr16_read.ATAC.bw
Cell1	H3K4me3	01	/storage/home/gzx103/scratch/S3V2norm_compare/hg38_EP/bw/c11_chr16_read.H3K4me3.bw
Cell2	ATAC	02	/storage/home/gzx103/scratch/S3V2norm_compare/hg38_cCREs/bw/c11_chr16_read.ATAC.bw
Cell2	H3K4me3	02	/storage/home/gzx103/scratch/S3V2norm_compare/hg38_EP/bw/c11_chr16_read.H3K4me3.bw
```

#####################################################################################

## How to run S3V2_IDEAS_ESMP pipeline
#### Use 'run_S3V2_IDEAS_ESMP.sh' to run S3norm pipeline.
##### After perparing the input data, user just need to set the parameters in "run_S3V2_IDEAS_ESMP.sh" to run S3V2_IDEAS_ESMP.
#####
```
### required inputs
###### the absolute path to the bin folder
script_dir='/storage/home/gzx103/group/software/S3V2_IDEAS_ESMP/bin/'
###### your output folder
output_dir='/storage/home/gzx103/scratch/S3V2norm_compare/S3V2_cCRE_pipeline_test/'
###### the absolute path to the your modified "metadata.for_master_peak_calls.txt" file
metadata='/storage/home/gzx103/scratch/S3V2norm_compare/S3V2_cCRE_pipeline_test/input_files/metadata.for_master_peak_calls.txt'
###### The output name
id_name='test_S3V2_IDEAS_ESMP_pipeline'


###### genome
GENOME='hg38'
###### genome size (can be found in the "S3V2_IDEAS_ESMP/genomesize/" folder)
GENOMESIZES='/storage/home/gzx103/group/software/S3V2_IDEAS_ESMP/genomesize/hg38.chrom.chr16.fortest.sizes'
###### blacklist (can be found in the "S3V2_IDEAS_ESMP/blacklist/" folder)
BLACK='/storage/home/gzx103/group/software/S3V2_IDEAS_ESMP/blacklist/hg38-blacklist.v2.bed'

###### number of threads in system
threads=4
###### bin size of the signal resolution
bin_size=200
###### email address
email='your_email@xxx.edu'


###### This one needs to be set to "T" for the first time of running the pipeline. 
get_sigtrack='T'
###### This one needs to be set to "T" for the first time of running the pipeline. 
normalization='T'
###### if user only want to get the S3V2 normalized bigwig files, set this to "T"
get_bw='T'
###### if user only want to run S3norm version2 without downstream epigenetic state call of hypersensitive state call, then this can be set to "F"
run_ideas='T'

```

##### Then Run:
```
time bash run_S3V2_IDEAS_ESMP.sh
```

#####################################################################################

## Outputs of S3V2_IDEAS_ESMP
### All outputs will be saved in the user provided "$output_dir".
##### The epigenetic state genome segmentation (multiple epigenetic features) and master peak list (one epigenetic feature) will be saved in the following folder: "your_output_name_IDEAS_output/"
```
ls -ltrh test_S3V2_IDEAS_ESMP_pipeline_IDEAS_output/
```
##### The S3norm normalized average read counts will be save in the "test_S3V2_IDEAS_cCRE_pipeline_bws_RC" folder
##### The -log10(p-value) based on S3norm normalized average read counts will be save in the "test_S3V2_IDEAS_cCRE_pipeline_bws_NBP" folder
##### The signal composition of the epigenetic state will be the "test_S3V2_IDEAS_cCRE_pipeline.pdf"
##### The genome segmentation will be saved in "test_S3V2_IDEAS_ESMP_pipeline_IDEAS_output/Tracks/" folder
##### If there is one epigenetic feature, a master peak list will be saved as the "test_S3V2_IDEAS_ESMP_pipeline.cCRE.M.bed" file

