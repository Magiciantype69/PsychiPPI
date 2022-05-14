import re
import sys
import os
import math
import random

print('ID geneset verify 1000000 times')
GeneID_match_Genename = {}
Human_cell_adhesion_molecules=[]; Actin_filament_bundle_assembly_gene=[]; Human_house_keeping_genes=[]; Neuronal_markers=[]; FMRP_target_in_mouse=[]; RBFOX_Splice_target=[]; G2C_human_PSD=[]; ARC_complex_genes=[];
NMDAR_complex_genes=[]; Constraint_genes=[]; Brain_preferentially_expressed=[]; Developmental_delay_genes=[]; ASD_genes =[]; Vulnerable_ASD_risk_genes=[]; ASD_DD_NEURO=[]; Neurodegenerative=[]
GenePool = []
with open(r'E:\OneDrive\学术\PPI article\data\geneset verify\geneset_verify_protein-coding.txt') as rd:
    text = rd.readline()
    text = text.strip()
    temp = text.split('\t')
    Human_cell_adhesion_molecules.append(temp[4]); Actin_filament_bundle_assembly_gene.append(temp[5]); Human_house_keeping_genes.append(temp[6]); Neuronal_markers.append(temp[12]); FMRP_target_in_mouse.append(temp[13]);
    RBFOX_Splice_target.append(temp[15]); G2C_human_PSD.append(temp[16]); ARC_complex_genes.append(temp[18]); NMDAR_complex_genes.append(temp[19]); Constraint_genes.append(temp[20]); Brain_preferentially_expressed.append(temp[21]);
    Developmental_delay_genes.append(temp[24]); ASD_genes.append(temp[26]); Vulnerable_ASD_risk_genes.append(temp[27]); ASD_DD_NEURO.append(temp[28]); Neurodegenerative.append(temp[29])
    for line in rd:
        line = line.strip()
        temp = line.split('\t')
        Genename = temp[1]
        GeneID = temp[0]
        GeneID_match_Genename[GeneID] = Genename
        GenePool.append(temp[0])
        if temp[4] =='1':
            Human_cell_adhesion_molecules.append(temp[0])
        if temp[5] =='1':
            Actin_filament_bundle_assembly_gene.append(temp[0])
        if temp[5] =='1':
            Human_house_keeping_genes.append(temp[0])
        if temp[12] =='1':
            Neuronal_markers.append(temp[0])
        if temp[13] =='1':
            FMRP_target_in_mouse.append(temp[0])
        if temp[15] =='1':
            RBFOX_Splice_target.append(temp[0])
        if temp[16] =='1':
            G2C_human_PSD.append(temp[0])
        if temp[18] =='1':
            ARC_complex_genes.append(temp[0])
        if temp[19] =='1':
            NMDAR_complex_genes.append(temp[0])
        if temp[20] =='1':
            Constraint_genes.append(temp[0])
        if temp[21] =='1':
            Brain_preferentially_expressed.append(temp[0])
        if temp[24] =='1':
            Developmental_delay_genes.append(temp[0])
        if temp[26] =='1':
            ASD_genes.append(temp[0])
        if temp[27] =='1':
            Vulnerable_ASD_risk_genes.append(temp[0])
        if temp[28] =='1':
            ASD_DD_NEURO.append(temp[0])
        if temp[29] =='1':
            Neurodegenerative.append(temp[0])
    verify_geneset = [Human_cell_adhesion_molecules,Actin_filament_bundle_assembly_gene,Human_house_keeping_genes,Neuronal_markers,FMRP_target_in_mouse,RBFOX_Splice_target,G2C_human_PSD,ARC_complex_genes,
    NMDAR_complex_genes,Constraint_genes,Brain_preferentially_expressed,Developmental_delay_genes,ASD_genes ,Vulnerable_ASD_risk_genes,ASD_DD_NEURO,Neurodegenerative]


random_geneset = []
Disease_geneIDset = []
Disease_geneIDset_length = 0
output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\8. Geneset verify\ID\heatmap.txt','w+')
output.write('Test_geneset'+'\t'+'Randon > Disease times'+'\t'+'In Geneset count'+'\t'+'Gene in set'+'\n')
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\8. Geneset verify\ID\ID_genelist.txt') as fd:
    title = fd.readline()
    for line in fd:
        line = line.strip()
        temp = line.split('\t')
        Disease_geneIDset.append(temp[1])
    Disease_geneIDset_length=len(Disease_geneIDset)

    while 1:
        random_geneset.append(random.sample(GenePool,Disease_geneIDset_length))
        if len(random_geneset)>1:
            break

    Disease_mapped_genelist = []
    Disease_mapped_countlist = []
    count_list = []

    for geneset in verify_geneset:
        inset = 0
        mapped_gene = []
        count_times = 0
        for selected_gene in Disease_geneIDset:
            if selected_gene in geneset:
                inset += 1
                mapped_gene.append(GeneID_match_Genename[selected_gene])
        for random_gene in random_geneset:
            random_inset = 0
            for ran_gene in random_gene:
                if ran_gene in geneset:
                    random_inset += 1
            if random_inset >= inset:
                count_times += 1
        output.write(geneset[0]+'\t'+str(count_times)+'\t'+str(inset)+'\t'+"\t".join(mapped_gene)+'\n')
    output.close()
