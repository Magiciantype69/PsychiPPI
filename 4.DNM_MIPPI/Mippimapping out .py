import re
import sys
import os
import math

print('birth defect out interface mutation is mapping mippi')
output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\4. DNM_MIPPI_mapping\birth defect_Mippi_INTERFACE_out.txt','w+')


with open(r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping\birth defect_DNM_Protein_mutation_mapping.txt') as root:
    root.readline()
    for line in root:
        line = line.strip()
        text = line.split('\t')
        mutation_info = text[6]; On_interface = int(text[12])

        if On_interface == 1:
            with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\mippi_filter.txt') as wood:
                wood.readline()
                for text in wood:
                    row = text.strip()
                    row = row.split('\t')
                    prediction = row[5];predict_score=row[6]; AAchange = row[2]
                    SourceID = row[0]; PartnerID = row[3]; Partnersymbol =row[4]; Sourcesymbol = row[1]
                    if AAchange in mutation_info:
                        output.write(text)

output.close()
