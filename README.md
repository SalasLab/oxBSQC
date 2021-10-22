# oxBSQC

## oxBScut

The function makes a recommendation for a 5hmC cut-off to eliminate unstable 5hmC probes on the Illumina HumanMethylation BeadChip and evaluate the information loss with the low-quality probes.

matrix_5hmc is a matrix or dataframe of 5hmC beta values.
poobha_5hmc is a matrix or dataframe of p-values with out-of-band array hybridization that is corresponding to the 5hmC beta matrix.


## oxBSlabel
Using internal normalization control probes (NORM_A/G/C/T) on Illumina HumanMethylation BeadChip to predict tandem bisulfite (BS) and oxidative bisulfite (oxBS) treated samples. 

Parameter "rgset" for the function should be a "RGChannelSet" with a pheno matrix with columns "Sample_ID" that is
unique to each sample and “Subject" that is unique to each paired subject.


## Installation

devtools::install_github("SalasLab/oxBSQC")


## Load library 
library(oxBSQC)


## Example
load("data/oxBScut_example.RDATA")

oxBScut(RCD_5hmC,RCD_poobha)

load("data/oxBSlabel_example.RDATA")

oxBSlabel(example_rgset)

