#!/bin/python
import os
import sys

wk="/work/long_lab/qli/hic_archieve/revision/LLR/GTEx/"
gtex_list = []
f=open(wk+"GTEx.list","r")
line=f.readline()
while line:
    gtex_list.append(line.strip())
    line=f.readline()
f.close()

sbatch_file=open(wk+"Permutation_sbatch.cmd","w")
for f_index in range(len(gtex_list)):
    if not os.path.exists(wk+gtex_list[f_index]+'/Permutation_1000'):
        cmd='mkdir -p '+wk+gtex_list[f_index]+'/Permutation_1000'
        os.system(cmd)
    
    cmd_file=open(wk+gtex_list[f_index]+"/Permutation_1000.cmd","w")
    cmd_file.write('#!/bin/bash\n')
    cmd_file.write('#SBATCH --job-name='+gtex_list[f_index]+'_Permutation_1000\n')
    cmd_file.write('#SBATCH --chdir='+wk+gtex_list[f_index]+'/Permutation_1000\n')
    cmd_file.write('#SBATCH --error=Permutation_1000.error\n')
    cmd_file.write('#SBATCH --out=Permutation_1000.out\n')
    cmd_file.write('#SBATCH --mem=10G\n')
    cmd_file.write('#SBATCH --nodes=1\n')
    cmd_file.write('#SBATCH --ntasks=1\n')
    cmd_file.write('#SBATCH --cpus-per-task=1\n')
    cmd_file.write('#SBATCH --time=7-00:00:00\n')
    cmd_file.write('#SBATCH --partition=single,lattice,parallel,theia,cpu2019,cpu2013\n')
    
    f_tag=open(wk+gtex_list[f_index]+".sig.genes.index.tsv","r")
    f_line=f_tag.readline()
    gene_written_list=[]
    while f_line:
        f_line_arr=f_line.split("\t")
        if f_line_arr[0] not in gene_written_list:
            cmd_file.write('java -jar /work/long_lab/qli/hic_archieve/revision/LLR/LLR.jar permutation -ipf '+wk+gtex_list[f_index]+"/"+gtex_list[f_index]+'.tsv.var.above.10.0.tsv.sorted -p 1000 -opf '+f_line_arr[0]+'_p1000 -index '+f_line_arr[2]+'\n')
            gene_written_list.append(f_line_arr[0])
        f_line=f_tag.readline()
    cmd_file.close()
    sbatch_file.write('sbatch '+wk+gtex_list[f_index]+"/Permutation_1000.cmd\n")
sbatch_file.close()

###MSSNG
#wk2='/work/long_lab/qli/hic_archieve/revision/LLR/MSSNG'
# 'java -jar /work/long_lab/qli/hic_archieve/revision/LLR/LLR.jar permutation -ipf '+wk2+'/pheno_0_1_adj.tsv -p 1000 -opf Permutation_1000/pheno_0_1_adj_p1000.tsv \n'