import re
import sys
import os
import math
from scipy import stats


output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\2. DNM_on_interface_enrichment\result\All enrichiment_result.txt','w+')
output.write('Disease'+'\t'+'PPI_mutation_count'+'\t'+'Mutation_count'+'\t'+'PPI_length_sum'+'\t'+'Protein_length_sum'+'\t'+ 'Domain_length_sum'+'\t'
            'interface_exact_binom_P_value'+'\t'+'interface enrichiment'+'\t'+'interface upper'+'\t'+'interface lower'+'\t'+
            'noninterface_exact_binom_P_value'+'\t'+'noninterface enrichiment'+'\t'+'noninterface upper'+'\t'+'noninterface lower'+
            '\t'+'domain_exact_binom_P_value'+'\t'+'domain enrichiment'+'\t'+'domain upper'+'\t'+'domain lower'+'\n')
rootdir = r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping'
name_list= os.listdir(rootdir)
for i in range(0,len(name_list)):
    path = os.path.join(rootdir,name_list[i])
    if os.path.isfile(path):
        Mutation_on_interface = 0
        Mutation_on_domain = 0
        count = 0
        All_Uniprotlist = []
        name = str(name_list[i])
        with open(path) as fp:
            fp.readline()
            for line in fp:
                line = line.strip()
                temp = line.split("\t")
                All_Uniprotlist.append(temp[0])
                Mutation_on_interface += int(temp[12])
                Mutation_on_domain += int(temp[13])
                count += int(temp[14])
        All_Uniprotlist = list(set(All_Uniprotlist))

        Uniprot_protein_length = {}
        Uniprot_PPI_length = {}
        Uniprot_domain_length = {}
        Uniprot_other_length = {}
        with open(r'E:\OneDrive\学术\PPI article\data\mapped data\2. DNM_on_interface_enrichment\Uniprot_Genename_length_domain_mapping.txt') as fp:
            fp.readline()
            for line in fp:
                line = line.strip()
                temp = line.split("\t")
                UniprotID = temp[0]
                Protein_length = temp[2]
                PPI_length = temp[3]
                Domain_length = temp[4]
                Other_length = temp[5]
                Uniprot_protein_length[UniprotID] = Protein_length
                Uniprot_PPI_length[UniprotID] = PPI_length
                Uniprot_domain_length[UniprotID] = Domain_length
                Uniprot_other_length[UniprotID] = Other_length


        All_protein_length = 0
        All_PPI_length = 0
        All_domain_length = 0
        for i in All_Uniprotlist:
            if i in Uniprot_protein_length:
                All_protein_length += int(Uniprot_protein_length[i])
                All_PPI_length += int(Uniprot_PPI_length[i])
                All_domain_length += int(Uniprot_domain_length[i])

        if Mutation_on_interface >0:
            M1 = Mutation_on_interface
            M2 = count - Mutation_on_interface - Mutation_on_domain
            M3 = Mutation_on_domain
            A_PPI = All_PPI_length
            A_domain = All_domain_length
            A_P = All_protein_length - All_PPI_length - All_domain_length
            A_length = All_protein_length

            interface_exact_binom = stats.binom_test(M1, count, A_PPI/A_length)
            noninterface_exact_binom = stats.binom_test(M2, count, (A_P/A_length))
            domain_exact_binom = stats.binom_test(M3, count, (A_domain/A_length))

            interface_enrichiment = (M1/count)/(A_PPI/A_length)
            noninterface_enrichiment = (M2/count)/(A_P/A_length)
            domain_enrichiment = (M3/count)/(A_domain/A_length)

            in_p1 =(1-(M1/count))/M1
            in_p2 = (1-(A_PPI/A_length))/A_PPI

            interface_upper = math.exp(math.log(interface_enrichiment)+1.96*(in_p1+in_p2)**0.5)
            interface_lower = math.exp(math.log(interface_enrichiment)-1.96*(in_p1+in_p2)**0.5)

            nonin_p1 =(1-(M2/count))/M2
            nonin_p2 = (1-(A_P/A_length))/A_P
            noninterface_upper = math.exp(math.log(noninterface_enrichiment)+1.96*(nonin_p1+nonin_p2)**0.5)
            noninterface_lower = math.exp(math.log(noninterface_enrichiment)-1.96*(nonin_p1+nonin_p2)**0.5)

            doin_p1 =(1-(M3/count))/M3
            doin_p2 = (1-(A_domain/A_length))/A_domain
            dointerface_upper = math.exp(math.log(domain_enrichiment)+1.96*(doin_p1+doin_p2)**0.5)
            dointerface_lower = math.exp(math.log(domain_enrichiment)-1.96*(doin_p1+doin_p2)**0.5)

            output.write(name+'\t'+str(Mutation_on_interface)+'\t'+str(count)+'\t'+str(All_PPI_length)+'\t'+str(All_protein_length)+'\t'+str(All_domain_length)+'\t'+
                        str(interface_exact_binom)+'\t'+str(interface_enrichiment)+'\t'+str(interface_upper)+'\t'+str(interface_lower)+'\t'+
                        str(noninterface_exact_binom)+'\t'+str(noninterface_enrichiment)+'\t'+str(noninterface_upper)+'\t'+str(noninterface_lower)+'\t'+
                        str(domain_exact_binom)+'\t'+str(domain_enrichiment)+'\t'+str(dointerface_upper)+'\t'+str(dointerface_lower)+'\t'+'\n')
output.close()
