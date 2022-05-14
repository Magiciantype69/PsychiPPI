import re
import sys
import os
import math

print('DNM_mapping_partner(ID) for degree')

output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\6. DNM_mapping_for_degree\SCZ.txt','w+')
output.write('P1'+'\t' + 'P2' + '\t'+ 'Source' +'\t'+ 'P1_symbol'+'\t' +'P2_symbol'+'\t'+'P1_entrez_id'+'\t'+'P2_entrez_id'+'\t'+'P1_PPI_mutation_num'+'\t'+'P2_PPI_mutation_num'+'\t'
+'mutation info'+'\n')

with open(r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping\psychiatric disorder\Schizophrenia (SCZ)_DNM_Protein_mutation_mapping.txt') as fp:
    fp.readline()
    for line in fp:

        Protein_Change_list = []
        line = line.strip()
        text = line.split('\t')
        UniprotID = text[0]; Protein_Change = text[11] ; Mutation_on_interface = int(text[12]) ; AAchange = text[6]
        if Mutation_on_interface == 1:
            Protein_change = Protein_Change.strip()
            Protein_Change_list = Protein_Change.split(',')
            Protein_Change_list = list(set(Protein_Change_list))


            with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\InterfacesHQ_with_detail.txt') as fd:
                fd.readline()
                for row in fd:
                    P1_PPI_mutation_num = 0
                    P2_PPI_mutation_num = 0
                    interfaces_P1 = []
                    interfaces_P2 = []
                    row = row.strip()
                    temp = row.split('\t')
                    P1 = temp[0]; P2 = temp[1] ; Source = temp[2] ; P1_symbol = temp[3] ;P2_symbol = temp[4] ; P1_entreid = temp[5] ; P2_entreid = temp[6] ;
                    P1_IRES = temp[7]; P2_IRES = temp[8] ; P1_protein_length = temp[9] ; P2_protein_length = temp[10] ; P1_binding_length = temp[11] ; P2_binding_length = temp[12]
                    if P1_IRES != '[]':
                        P1_IRES=P1_IRES[1:-1].split(',')
                        interfaces_P1=list(map(int,P1_IRES))
                    else:
                        interfaces_P1 = []

                    if P2_IRES != '[]':
                        P2_IRES=P2_IRES[1:-1].split(',')
                        interfaces_P2=list(map(int,P2_IRES))
                    else:
                        interfaces_P2 = []

                    if UniprotID == P1:
                        for i in Protein_Change_list:
                            mutation_site = re.findall(r'\d+',i)
                            if int(mutation_site[0]) in interfaces_P1:
                                P1_PPI_mutation_num = 1

                        if P1_PPI_mutation_num == 1:
                            output.write(P1+'\t' + P2 + '\t'+ Source +'\t'+ P1_symbol+'\t' +P2_symbol+'\t'+P1_entreid+'\t'+P2_entreid+'\t'+str(P1_PPI_mutation_num)+'\t'+str(P2_PPI_mutation_num)+'\t' + AAchange +'\n')
                    elif UniprotID == P2:
                        for i in Protein_Change_list:
                            mutation_site = re.findall(r'\d+',i)
                            if int(mutation_site[0]) in interfaces_P2:
                                P2_PPI_mutation_num = 1
                        if P2_PPI_mutation_num == 1:
                            output.write(P1+'\t' + P2 + '\t'+ Source +'\t'+ P1_symbol+'\t' +P2_symbol+'\t'+P1_entreid+'\t'+P2_entreid+'\t'+str(P1_PPI_mutation_num)+'\t'+str(P2_PPI_mutation_num)+'\t' + AAchange +'\n')
output.close()
