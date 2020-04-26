# 3DMM User Manual
3DMM stands for 3D Mixed Model which performs association study between genotype and phenotype using a transformed linear mixed model and integrating 3D genotype information.
Here, we present java code and executable jar file to run simulation study for 3DMM project. The association study is completed by another tool named JAWAMix5 developed in my lab. Details for JAMAMix5 can be found through following URL: https://github.com/theLongLab/Jawamix5.git.
Note: The method 3DMM is also known as "compound" method in JAMAMix5.
## Contents
1. [Installation](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#installation)
2. [Functions](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#functions)
3. [Commands](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#commands)
4. [Input/output file format](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#inputoutput-file-format)
5. [Acknowledgements](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#acknowledgements)

## Installation
We suggest to download the 3DMM_sim.jar file directly and run it using the regular command for jar file:
  java -jar /path/to/3DMM_sim.jar
Otherwise, please clone the source code to your java development environment such as Eclipse.
## Functions
We provide two functions:
1. simulation: to generate simulated phenotypes
2. power: to calculate power of the association method based on its p-values.

## Commands
1. simulation
  * -genotype: input genotype file in HDF5 format
  * -ov: output file recording causal variants
  * -oc: output file recording causal paired-regions
  * -op: output file for phenotypes
  * -sr: simulation phenotype report
  * -info: position information for paired-regions
  * -type: compound/SNP
  * -nr: a numerical value for number of rounds
  * -nc: a numerical value for number of causal paired-regions
  * -ns: a numerical value for number of causal variants in each part of paired-regions
  * -f: a numerical value for length of flank region
  * -min: a numerical value for minimal MAF of SNPs to be used
  * -max: a numerical value for maximal MAF of SNPs to be used
  * -gc: a numerical value for heritability of causal paired-regions used
  * -p: cross-region pattern (add/hetero/both/compensate)
  * -th: a numerical value indicating liability threshold (optional)
  * -binary: a Boolean value indicate binary phenotype (false) and quantitative phenotype (true)
2. power
  * -p_path_p: the file contains all p-values of a certain model
  * -padj: the file containing the cut-off for significant p-values
  * -model: compound/local/emmax/skat
  * -used_var: the file containing all causal variants
  * -of: output folder
  * -cre: criteria: 0: for compound method, 1: one part of the paired-regions, 2: two parts of the paired-regions
  * -w: length of a window in genome

## Input/output file format
**genotype** and **phenotype** please refer to those used in JAWAMix5

**position information for paired-regions** is in bellow format
  #Index  part1_chr;part1_start;part1_end part2_chr;part2_start;part2_end
  C0      1;840000;850000 1;890000;900000

**p_path_p** contains all path for p-values generating from each round
  /path/to/p-values/of/1st/round
  /path/to/p-values/of/2nd/round
  e.g.
  ./Compound_VO.0.Phenotype_round_0.csv
  ./Compound_VO.1.Phenotype_round_1.csv 

**padj**  has the cut-off for significant paired-regions
  compound,2.124E-6
  local,5.91422E-8
  skat,4.167E-8
  emmax,4.26E-10
## Acknowledgements
I appreciate to have such an opportunity to leading this 3DMM project. I am grateful for all supports provided by my supervisor Dr. Long and my lab mates.
