import re
import sys
import os
import math
from pyfaidx import Fasta


#mapping UniprotID with each Protein interface region
UniprotID_interface = {}
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\Uniprot_Genename_All_PPI_sites.txt') as fp:
    fp.readline()
    for line in fp:
        line= line.strip('\n')
        temp = line.split("\t")
        if len(temp)<=3:
            UniprotID = temp[0]
        if len(temp) <3:
            interface_residue = ''
        else:
            interface_residue = temp[2]
        UniprotID_interface[UniprotID] = interface_residue

#mapping UniprotID with each Protein domain region
UniprotID_domain = {}
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\Uniprot_Genename_All_domain_sites.txt') as fx:
    fx.readline()
    for line in fx:
        line= line.strip('\n')
        temp = line.split("\t")
        if len(temp)<=3:
            UniprotID = temp[0]
        if len(temp) <3:
            domain_residue = ''
        else:
            domain_residue = temp[2]
        UniprotID_domain[UniprotID] = domain_residue


#extract all crucial information from annotation masterfile
Disease_category = {}
Disease_list = []
Category_list =[]
DNM_PPI_interface_mapping = open(r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping\All_DNM_Protein_mutation_mapping.txt','w+')
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\psymukb filter.txt') as fp:
    text = fp.readline()
    text = text.strip('\n')
    title_text = text+'\t'+'Mutation_on_interface'+'\t'+'Mutation_on_domain'+'\t'+'count'+'\n'
    DNM_PPI_interface_mapping.write(title_text)
    for line in fp:
        line = line.strip('\n')
        temp = line.split("\t")
        UniprotID = temp[0]; GeneID = temp[1];Gene_refGene = temp[2] ; ExonicFunc_refGene = temp[3] ; PrimaryPhenotype = temp[4] ; DisorderCategory = temp[5] ; AAChange_refGene =temp[6]; SIFT_score = temp[7] ; Polyphen2_HDIV_score = temp[8] ; CADD_phred = temp[9];
        Protein_Change = temp[11];
        Disease_category[PrimaryPhenotype] = DisorderCategory
        Disease_list.append(PrimaryPhenotype)
        Category_list.append(DisorderCategory)

        Protein_Change_list = Protein_Change.split(',')
        Protein_Change_residue = []
        pos = []

        for residue in Protein_Change_list:
            if residue not in Protein_Change_residue:
                Protein_Change_residue.append(residue)
                pos_temp = re.findall('\d+',residue)
                pos.append(pos_temp[0])

        Mutation_on_interface = 0
        Mutation_on_domain = 0
        if UniprotID in UniprotID_interface:
            interface_list = UniprotID_interface[UniprotID].split(',')
            interface_list=[x.strip() for x in interface_list if x.strip() != '']
            residue_list= list(map(int,interface_list))
            if len(residue_list) <1:
                residue_list = [0]
            for Protein in pos:
                if int(Protein) in residue_list:
                    Mutation_on_interface = 1
            if Mutation_on_interface == 0:
                if UniprotID in UniprotID_domain:
                    domain_list = UniprotID_domain[UniprotID].split(',')
                    domain_list=[x.strip() for x in domain_list if x.strip() != '']
                    domain_residue_list= list(map(int,domain_list))
                    if len(domain_residue_list) <1:
                        domain_residue_list = [0]
                    for Protein in pos:
                        if int(Protein) in domain_residue_list:
                            Mutation_on_domain =1


        DNM_PPI_interface_mapping.write(line +'\t' +str(Mutation_on_interface)+'\t' +str(Mutation_on_domain)+'\t' +'1'+'\n')
    Disease_list = list(set(Disease_list))
    Category_list = list(set(Category_list))
DNM_PPI_interface_mapping.close()

#split the information by Disease_list
for disease in Disease_list:
    dir1 = r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping'+'\\'
    dir2 = r'_DNM_Protein_mutation_mapping.txt'
    disease_name = disease.replace('/','')
    dir = dir1+Disease_category[disease]+'\\'+disease_name+dir2
    disease_output = open(dir,'w+')
    disease_output.write(title_text)
    with open(r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping\All_DNM_Protein_mutation_mapping.txt') as fd:
        text = fd.readline()
        for line in fd:
            line = line.strip()
            temp = line.split("\t")
            if temp[4] == disease:
                disease_output.write(line+'\n')

    disease_output.close()

for category in Category_list:
    dir1 = r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping'+'\\'
    dir2 = r'_DNM_Protein_mutation_mapping.txt'
    dir = dir1+ category + dir2
    category_output = open(dir,'w+')
    category_output.write(title_text)
    with open(r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping\All_DNM_Protein_mutation_mapping.txt') as fd:
        text = fd.readline()
        for line in fd:
            line = line.strip()
            temp = line.split("\t")
            if temp[5] == category:
                category_output.write(line+'\n')

    category_output.close()
