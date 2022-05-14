import re
import sys
import os
import math

print('ID_only.txt interface mutation is mapping mippi')
output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\4. DNM_MIPPI_mapping\Schizophrenia_Mippi_INTERFACE.txt','w+')

mapped = 0
out = 0
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\6. DNM_mapping_for_degree\SCZ_only.txt') as root:
    output.write((root.readline()).strip('\n')+'\t'+'prediction'+'\t'+'predict score'+'\n')
    for line in root:
        line = line.strip()
        text = line.split('\t')
        P1_symbol=text[3];P2_symbol=text[4];P1_entrez_id=text[5];P2_entrez_id=text[6];
        P1_PPI_mutation_num=text[7];P2_PPI_mutation_num=text[8];mutation_info=text[9];
        found = 0

        with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\mippi_filter.txt') as wood:
            wood.readline()
            for row in wood:
                row = row.strip()
                row = row.split('\t')
                prediction = row[5];predict_score=row[6]; AAchange = row[2]
                SourceID = row[0]; PartnerID = row[3]; Partnersymbol =row[4]; Sourcesymbol = row[1]
                if AAchange in mutation_info:
                    if P1_PPI_mutation_num == '1':
                        mutation_geneID = P1_entrez_id
                        mutation_partnerID = P2_entrez_id
                    elif P2_PPI_mutation_num =='1':
                        mutation_geneID = P2_entrez_id
                        mutation_partnerID = P1_entrez_id

                    if mutation_geneID == SourceID and mutation_partnerID == PartnerID:
                        output.write(line+'\t'+prediction+'\t'+predict_score+'\n')
                        found = 1
                        mapped += 1
                        break
        if found  == 0:
            output.write(line+'\t'+'out of prediction'+'\n')
            out += 1
output.close()

print('mapped PPI :' +str(mapped))
print('out of mapping:'+str(out))
