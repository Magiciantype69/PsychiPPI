import re
import sys
import os
import math

Uniprot_domainlength = {}

with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\Uniprot_Genename_All_domain_sites.txt') as fp:
    fp.readline()
    for line in fp:
        line=line.strip('\n')
        text=line.split('\t')
        UniprotID = text[0]
        if text[2]=='':
            Uniprot_domainlength[UniprotID]=0
        else:
            domainlength = len(text[2].split(','))
            Uniprot_domainlength[UniprotID] = domainlength

output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\2. DNM_on_interface_enrichment\Uniprot_Genename_length_domain_mapping.txt','w+')
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\2. DNM_on_interface_enrichment\Uniprot_Genename_length_mapping.txt') as fx:
    output.write('UniprotID'+'\t'+'Genename'+'\t'+'Protein_length'+'\t'+'PPI length'+'\t'+'Domain_length'+'\t'+'Other_length'+'\n')
    fx.readline()
    for line in fx:
        line = line.strip('\n')
        text = line.split('\t')
        Uniprot= text[0]; Genename = text[1]
        Protein_length =int(text[2])
        PPI_length = int(text[3])
        Domain_length = int(Uniprot_domainlength[Uniprot])
        Other_length = Protein_length - PPI_length -Domain_length
        output.write(Uniprot+'\t'+Genename+'\t'+str(Protein_length)+'\t'+str(PPI_length)+'\t'+str(Domain_length)+'\t'+str(Other_length)+'\n')
output.close()
