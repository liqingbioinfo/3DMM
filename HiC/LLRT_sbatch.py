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

sbatch_file=open(wk+"LLRT_sbatch.cmd","w")

for f_index in range(len(gtex_list)):
    if not os.path.exists(wk+gtex_list[f_index]+'/LLR_br'):
        cmd='mkdir -p '+wk+gtex_list[f_index]+'/LLR_br'
        os.system(cmd)
    if not os.path.exists(wk+gtex_list[f_index]+'/LLR_sr'):
        cmd='mkdir -p '+wk+gtex_list[f_index]+'/LLR_sr'
        os.system(cmd)
    
    tissue_sbatch_file_br=open(wk+gtex_list[f_index]+"/"+gtex_list[f_index]+"_LLR_sbatch_br.cmd","w")
    tissue_sbatch_file_sr=open(wk+gtex_list[f_index]+"/"+gtex_list[f_index]+"_LLR_sbatch_sr.cmd","w")
    f_tag=open(wk+gtex_list[f_index]+".sig.genes.index.tsv","r")
    f_line=f_tag.readline()
    gene_analyzed_list=[]
    while f_line:
        f_line_arr=f_line.split("\t") # gene_name compound_name gene_index
        gene_name=f_line_arr[0]
        if gene_name not in gene_analyzed_list:
            if not os.path.exists(wk+gtex_list[f_index]+'/LLR_br/'+gene_name):
                cmd='mkdir -p '+wk+gtex_list[f_index]+'/LLR_br/'+gene_name
                os.system(cmd)
            
            cmd_file=open(wk+gtex_list[f_index]+'/LLR_br/'+gene_name+'.br.cmd','w')
            cmd_file.write('#!/bin/bash\n')
            cmd_file.write('#SBATCH --job-name='+gene_name+'_br_LLR\n')
            cmd_file.write('#SBATCH --chdir='+wk+gtex_list[f_index]+'/LLR_br/\n')
            cmd_file.write('#SBATCH --error='+gene_name+'_br_LLR.error\n')
            cmd_file.write('#SBATCH --out='+gene_name+'_br_LLR.out\n')
            cmd_file.write('#SBATCH --mem=10G\n')
            cmd_file.write('#SBATCH --nodes=1\n')
            cmd_file.write('#SBATCH --ntasks=1\n')
            cmd_file.write('#SBATCH --cpus-per-task=1\n')
            cmd_file.write('#SBATCH --time=7-00:00:00\n')
            cmd_file.write('#SBATCH --partition=single,lattice,parallel,theia,cpu2019,cpu2013\n')
            cmd_file.write('java -jar /work/long_lab/qli/hic_archieve/revision/LLR/LLR.jar LLRT -ig '+wk+gtex_list[f_index]+"/"+gtex_list[f_index]+'.hdf5 -ip '+wk+gtex_list[f_index]+'/Permutation_1000/'+gene_name+'_p1000 -ik_g '+wk+gtex_list[f_index]+"/"+gtex_list[f_index]+'.rescaled.IBS -of ./'+gene_name+'/ -ic '+wk+gtex_list[f_index]+'/GTEx_sig_res_'+gtex_list[f_index]+'.txt.br.c\n')
            cmd_file.close()
            tissue_sbatch_file_br.write('sbatch '+wk+gtex_list[f_index]+'/LLR_br/'+gene_name+'.br.cmd\n')
            
            if not os.path.exists(wk+gtex_list[f_index]+'/LLR_sr/'+gene_name):
                cmd='mkdir -p '+wk+gtex_list[f_index]+'/LLR_sr/'+gene_name
                os.system(cmd)
            
            cmd_file=open(wk+gtex_list[f_index]+'/LLR_sr/'+gene_name+'.sr.cmd','w')
            cmd_file.write('#!/bin/bash\n')
            cmd_file.write('#SBATCH --job-name='+gene_name+'_sr_LLR\n')
            cmd_file.write('#SBATCH --chdir='+wk+gtex_list[f_index]+'/LLR_sr/\n')
            cmd_file.write('#SBATCH --error='+gene_name+'_sr_LLR.error\n')
            cmd_file.write('#SBATCH --out='+gene_name+'_sr_LLR.out\n')
            cmd_file.write('#SBATCH --mem=10G\n')
            cmd_file.write('#SBATCH --nodes=1\n')
            cmd_file.write('#SBATCH --ntasks=1\n')
            cmd_file.write('#SBATCH --cpus-per-task=1\n')
            cmd_file.write('#SBATCH --time=7-00:00:00\n')
            cmd_file.write('#SBATCH --partition=single,lattice,parallel,theia,cpu2019,cpu2013\n')  
            cmd_file.write('java -jar /work/long_lab/qli/hic_archieve/revision/LLR/LLR.jar LLRT -ig '+wk+gtex_list[f_index]+"/"+gtex_list[f_index]+'.hdf5 -ip '+wk+gtex_list[f_index]+'/Permutation_1000/'+gene_name+'_p1000 -ik_g '+wk+gtex_list[f_index]+"/"+gtex_list[f_index]+'.rescaled.IBS -of ./'+gene_name+'/ -ic '+wk+gtex_list[f_index]+'/GTEx_sig_res_'+gtex_list[f_index]+'.txt.sr.c\n')
            cmd_file.close()
            tissue_sbatch_file_sr.write('sbatch '+wk+gtex_list[f_index]+'/LLR_sr/'+gene_name+'.sr.cmd\n')
            
            gene_analyzed_list.append(gene_name)
        f_line=f_tag.readline()
        
    tissue_sbatch_file_br.close()
    tissue_sbatch_file_sr.close()
    f_tag.close()
    sbatch_file.write('bash '+wk+gtex_list[f_index]+"/"+gtex_list[f_index]+"_LLR_sbatch_br.cmd\n")
    sbatch_file.write('bash '+wk+gtex_list[f_index]+"/"+gtex_list[f_index]+"_LLR_sbatch_sr.cmd\n")
sbatch_file.close()
            