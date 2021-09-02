# oxBSlabel

Using internal normalization control probes (NORM_A/G/C/T) on Illumina HumanMethylation BeadChip to predict tandem bisulfite (BS) and oxidative bisulfite (oxBS) treated samples. 


## Installation

devtools::install_github("SalasLab/oxBSlabel")


## Load library 
library(oxBSlabel)


## Example
load("data/example_rgset.RDATA")
oxBSlabel(example_rgset)



