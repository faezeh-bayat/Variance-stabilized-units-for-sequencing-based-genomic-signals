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
install.packages('bigmemory')
install.packages('pracma')

bedtools
######(https://bedtools.readthedocs.io/en/latest/content/installation.html)
```

## Installing VSS pipeline
```
git clone https://github.com/faezeh-bayat/Variance-stabilized-units-for-sequencing-based-genomic-signals.git
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

## Output of VSS pipeline
#### VSS generates variance-stabilized signals. The output of the pipeline is a bedGraph format transformed signal.

## How to run VSS pipeline
#### 1. Provide inputs in VSS.sh
```
### path to the input files
input_directory="/absolute_path_of_the_users_input_files"
### Indicate whether you want to build your own model("user_specified") or use the default model ("default")  
variance_stabilization_model="user_specified" or "default"
### Two input replicates for building the mean-varaince relationship model
replicate1_signals_for_training_the_model="rep1.bedGraph" or "False" if the variance_stabilization_model is "default"
replicate2_signals_for_training_the_model="rep2.bedGraph" or "False" if the variance_stabilization_model is "default"
### Indicate which chromosomes you want to build the model on 
chromosomes_to_build_the_model="chr21" or "False" if the variance_stabilization_model is "default". Chromosomes can also be in "all" format if user wants to build model for all chromosomes.
###
signals_to_be_variance_stabilized="rep1.bedGraph"
### Indicate which chromosome of untransformed signals you want to transform
chromosomes_to_be_stabilized="chr21" ###Chromosomes can also be in "all" format if user wants to get all chromosomes' variance-stabilized signals
###
Source_directory="/absolute_path_where_VSS_is_installed"

module load r
module load r/3.4.4
module load bedtools/2.27.1

cd "$Source_directory"
cd "Source"
Rscript VSS.R $input_directory \
                $variance_stabilization_model \
                $replicate1_signals_for_training_the_model \
                $replicate2_signals_for_training_the_model \
                $chromosomes_to_build_the_model \
                $signals_to_be_variance_stabilized \
                $chromosomes_to_be_stabilized \
                $Source_directory

```
#### 2. Run VSS.sh

```
chmode +x VSS.sh
./VSS.sh

```
#### 3. Variance-stabilized signals will be saved in the output folder at user provided "$input_directory".


##################################################################################################

#### You can accesee the "Variance-stabilized units for sequencing-based genomic signals" manuscript in:
https://www.biorxiv.org/content/10.1101/2020.01.31.929174v2

