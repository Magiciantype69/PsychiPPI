import re
import sys
import os
import math
from pyfaidx import Fasta

#extract all crucial information from annotation masterfile
DNM_PPI_interface_mapping = open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\psymukb filter1.txt','w+')
with open(r'E:\OneDrive\学术\PPI article\data\de novo mutation data\20190514_MasterFile_allDNMs-codingLoc_v1.5.txt') as fp:
    text = fp.readline()
    text = text.strip()
    title = text.split('\t')
    title_text = title[0] + '\t' + title[8] + '\t' + title[10] + '\t' + title[17] + '\t' + title[18] + '\t' + title[23] + '\t' + title[57] + '\t' + title[60] +'\t' +title[94]+'\t'+title[152]+'\t'+'Protein Change' +'\n'
    DNM_PPI_interface_mapping.write(title_text)
    for line in fp:
        Protein_Change = ''
        line = line.strip()
        temp = line.split("\t")
        GeneID = temp[0];Gene_refGene = temp[8] ; ExonicFunc_refGene = temp[10] ; PrimaryPhenotype = temp[17] ; DisorderCategory = temp[18] ; SIFT_score = temp[57] ; Polyphen2_HDIV_score = temp[60] ; CADD_phred = temp[94]
        AAChange_refGene = temp[23]; Haploinsufficiency_Score=temp[152]
        if SIFT_score == '':
            SIFT_score ='.'
        if Polyphen2_HDIV_score == '':
            Polyphen2_HDIV_score ='.'
        if CADD_phred == '':
            CADD_phred ='.'
        if Haploinsufficiency_Score == '':
            Haploinsufficiency_Score= '.'

        AAchange = AAChange_refGene.split(',' and ':')

        for element in AAchange:
            if 'p.' in element:
                if ',' in element:
                    element = element.split(',')
                    element = element[0]
                Protein_Change = Protein_Change +','+ element
                if Protein_Change.endswith('"'):
                    Protein_Change = Protein_Change[:-1]
                if Protein_Change[0] ==',':
                    Protein_Change = Protein_Change[1:]

        Protein_Change_list = Protein_Change.split(',')
        Protein_Change_residue = []
        pos = []
        if ExonicFunc_refGene == 'nonsynonymous SNV':
            for residue in Protein_Change_list:
                if residue not in Protein_Change_residue:
                    Protein_Change_residue.append(residue)
                    pos_temp = re.findall('\d+',residue)
                    pos.append(pos_temp[0])

            if ';' in Gene_refGene:
                if ';' in Gene_refGene:
                    Gene_refGene = Gene_refGene.split(';')
                    GeneID = GeneID.split(';')
                DNM_PPI_interface_mapping.write(GeneID[0] + '\t' + Gene_refGene[0] + '\t' + ExonicFunc_refGene + '\t' + PrimaryPhenotype + '\t' + DisorderCategory + '\t' + AAChange_refGene + '\t' + SIFT_score + '\t' +
                Polyphen2_HDIV_score + '\t' + CADD_phred+'\t'+ Haploinsufficiency_Score+'\t'+ Protein_Change +'\n')
                Gene_refGene = Gene_refGene[1]
                GeneID = GeneID[1]


            DNM_PPI_interface_mapping.write(GeneID + '\t' + Gene_refGene + '\t' + ExonicFunc_refGene + '\t' + PrimaryPhenotype + '\t' + DisorderCategory + '\t' + AAChange_refGene + '\t' + SIFT_score + '\t' +
            Polyphen2_HDIV_score + '\t' + CADD_phred+'\t'+ Haploinsufficiency_Score+'\t'+ Protein_Change +'\n')
DNM_PPI_interface_mapping.close()

Uniprot_GeneID ={}
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\psymukb filter old.txt') as fe:
    fe.readline()
    for line in fe:
        line = line.strip()
        line = line.split('\t')
        Uniprot = line[0]
        GeneID = line[1]
        Uniprot_GeneID[GeneID] = Uniprot

output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\psymukb filter.txt','w+')
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\psymukb filter1.txt') as fa:
    output.write('Uniprot'+'\t'+fa.readline())
    for line in fa:
        text = line.strip()
        text = text.split('\t')
        geneid = text[0]
        if geneid in Uniprot_GeneID:
            output.write(Uniprot_GeneID[geneid]+'\t'+line)
output.close()
