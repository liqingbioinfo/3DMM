# ILMM User Manual
ILMM stands for Interaction-Integrated Linear Mixed Model which performs association study between genotype and phenotype using a transformed linear mixed model and integrating 3D genotype information.
Here, we present java source code and executable jar file to run simulation study (ILMM_sim.jar) and likelihood ratio test (LLRT.jar) for ILMM project. The association study is completed by another tool named JAWAMix5 developed in my supervisor's lab. Details for JAMAMix5 can be found through following URL: https://github.com/theLongLab/Jawamix5.git.
Note: The method ILMM is also known as "compound" method in JAMAMix5.
## Contents
1. [Installation](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#installation)
2. [Functions](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#functions)
3. [Commands](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#commands)
4. [Input/output file format](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#inputoutput-file-format)
5. [Examples](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#examples)
6. [Acknowledgements](https://github.com/liqingbioinfo/3DMM/blob/master/README.md#acknowledgements)

## Installation
We suggest to download the ILMM_sim.jar and LLRT.jar file directly and run it using the regular command for jar file:
  java -jar /path/to/ILMM_sim.jar
  java -jar /path/to/LLRT.jar
Otherwise, please clone the source code to your java development environment such as Eclipse.

## Functions
We provide two functions:
### Simulation (ILMM_sim.jar)
1. simulation: to generate simulated phenotypes
2. power: to calculate power of the association method based on its p-values.
### LLRT (LLRT.jar)
1. permutation: generate distribution under null hypothesis
2. LLRT: likelihood ratio test
3. LLRT_p: get pvalue for likelihood ratio tests

## Commands
### Simulation (ILMM_sim.jar)
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
### LLRT (LLRT.jar)
1. permutation
  * -ipf: phenotype file name (seperated by \t)
  * -p: rounds of permutation (>100 suggested)
  * -opf: output phenotype file name (seperated by \t)
  * -index: the index of columns of phenotype to be permuated (df=0, useful while original phenotype file contains multiple columns) [optional]
2. LLRT
  * -ig: genotype file (hdf5)
  * -ip: phenotype file name (one file contains permuated phenotypes)
  * -ik_g: global kinship matrix
  * -of: output folder name
  * -ic: file name which contains compounds
  * -index: the index of columns of phenotype to be analyzed (df=0, useful while original phenotype file contains multiple columns) [optional]
3. LLRT_p
  * -io_br: input folder path (all results for one phenoytpe with double regions)
  * -io_sr: input folder path (all results for one phenoytpe with single region)
  * -ic: one compound name
  * -pheno: phenotpe file name
  * -p: rounds of permutation
  * -op: file name for output pvalue

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

## Examples
### Simulation (ILMM_sim.jar)
1. simulation
java -jar /path/to/ILMM_sim.jar simulation -genotype example.hdf5 -ov ./used_vars.txt -oc ./used_compounds.txt -op ./sim_pheno -sr ./sim_repot -info /path/to/compound_file -type compound -nr 500 -nc 1 -ns 5 -f 25000 -min 0.15 -max 0.49 -gc 0.3 -p add
2. power
java -jar /path/to/ILMM_sim.jar power -p_path_p compound_p_files.sorted.txt -padj ../p_adjust_value.txt -model compound -used_var used_vars.txt.balance.vars.txt -of ./power2/ -cre 0 -w 0
### LLRT (LLRT.jar)
1. permutation
java -jar /path/to/LLR.jar permutation -ipf example.tsv -p 1000 -opf example_p1000.tsv -index 0
2. LLRT
java -jar /path/to/LLR.jar LLRT -ig example.hdf5 -ip example_p1000.tsv -ik_g example.rescaled.IBS -of /path/to/br/ -ic example.txt.br.c
java -jar /path/to/LLR.jar LLRT -ig example.hdf5 -ip example_p1000.tsv -ik_g example.rescaled.IBS -of /path/to/sr/ -ic example.txt.sr.c
3. LLRT_p
java -jar /path/to/LLR.jar LLRT_p -io_br /path/to/br/ -io_sr /path/to/sr/ -ic 1:2350000:2360000_1:2480000:2490000 -pheno pheno_name -p 1000 -op_file pvalue
4. Note: To run multiple phenotypes in HPC, please adapt Permutation_sbatch.py, LLRT_sbatch.py, LLRT_p_sbatch.py to your own work path.

## Acknowledgements
I appreciate to have such a great opportunity to lead this ILMM project. I am grateful for all supports provided by my supervisor Dr. Long and my lab mates.
