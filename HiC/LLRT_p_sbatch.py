#!/bin/python
import os
import sys

##For tissue 1-10, p1000 together for sr and br
wk="/work/long_lab/qli/hic_archieve/revision/LLR/GTEx/"
gtex_list = []
f=open(wk+"GTEx.list","r")
line=f.readline()
while line:
    gtex_list.append(line.strip())
    line=f.readline()
f.close()

sbatch_file=open(wk+"LLRT_p_sbatch.cmd","w")

for f_index in range(len(gtex_list)):
    if not os.path.exists(wk+gtex_list[f_index]+'/LLR_pvalue'):
        cmd='mkdir -p '+wk+gtex_list[f_index]+'/LLR_pvalue'
        os.system(cmd)
    
    cmd_file=open(wk+gtex_list[f_index]+'/LLR_pvalue/LLRT_pvalue.cmd','w')
    cmd_file.write('#!/bin/bash\n')
    cmd_file.write('#SBATCH --job-name='+gtex_list[f_index]+'_LLRT_pvalue\n')
    cmd_file.write('#SBATCH --chdir='+wk+gtex_list[f_index]+'/LLR_pvalue\n')
    cmd_file.write('#SBATCH --error='+gtex_list[f_index]+'_LLRT_pvalue.error\n')
    cmd_file.write('#SBATCH --out='+gtex_list[f_index]+'_LLRT_pvalue.out\n')
    cmd_file.write('#SBATCH --mem=10G\n')
    cmd_file.write('#SBATCH --nodes=1\n')
    cmd_file.write('#SBATCH --ntasks=1\n')
    cmd_file.write('#SBATCH --cpus-per-task=1\n')
    cmd_file.write('#SBATCH --time=7-00:00:00\n')
    cmd_file.write('#SBATCH --partition=single,lattice,parallel,theia,cpu2019,cpu2013\n')
    
    f_tag=open(wk+gtex_list[f_index]+".sig.genes.index.tsv","r")
    f_line=f_tag.readline()
    while f_line:
        f_line_arr=f_line.split("\t") # gene_name compound_name gene_index
        cmd_file.write('java -jar /work/long_lab/qli/hic_archieve/revision/LLR/LLR.jar LLRT_p -io_br '+wk+gtex_list[f_index]+'/LLR_br/'+f_line_arr[0]+'/ -io_sr '+wk+gtex_list[f_index]+'/LLR_sr/'+f_line_arr[0]+'/ -ic '+f_line_arr[1]+' -pheno '+f_line_arr[0]+' -p 1000 -op_file '+f_line_arr[0]+'_'+f_line_arr[1]+'\n')
        f_line=f_tag.readline()
    f_tag.close()
    cmd_file.close()
    sbatch_file.write('sbatch '+wk+gtex_list[f_index]+'/LLR_pvalue/LLRT_pvalue.cmd\n')
    
sbatch_file.close()
        