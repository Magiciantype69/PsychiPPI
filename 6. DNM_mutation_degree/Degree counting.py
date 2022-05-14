import re
import sys
import os
import math

print('DNM_mapping_partner(ID) for degree counting')

output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\6. DNM_mapping_for_degree\Degree counting\ID counting degree.txt','w+')

output.write('UniprotID'+'\t'+'Mutation'+'\t'+'Degree'+'\n')

with open(r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping\psychiatric disorder\Intellectual disability (ID)_DNM_Protein_mutation_mapping.txt') as fp:
    fp.readline()
    for line in fp:
        line = line.strip()
        text = line.split('\t')
        UniprotID = text[0]; Protein_Change = text[11] ; Mutation_on_interface = int(text[12]) ; AAchange = text[6]
        Degree = 0
        if Mutation_on_interface == 1:
            with open(r'E:\OneDrive\学术\PPI article\data\mapped data\6. DNM_mapping_for_degree\ID_only.txt') as root:
                root.readline()
                for row in root:
                    row = row.strip()
                    temp = row.split('\t')
                    mutation_info = temp[9]
                    if AAchange == mutation_info:
                        Degree += 1
            if Degree != 0:
                output.write(UniprotID + '\t' + AAchange + '\t' + str(Degree)+'\n')

output.close
