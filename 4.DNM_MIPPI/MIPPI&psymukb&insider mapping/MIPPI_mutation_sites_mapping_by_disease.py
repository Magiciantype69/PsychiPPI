import re
import sys
import os
import math
from scipy import stats

print("MIPPI_mutation_sites_mapping_by_disease")
Disease_list = []
Category_list =[]
Disease_category = {}
with open(r'E:\OneDrive\学术\PPI article\data\de novo mutation data\20190514_MasterFile_allDNMs-codingLoc_v1.5.txt') as fp:
    text = fp.readline()
    for line in fp:
        line = line.strip()
        temp = line.split("\t")
        PrimaryPhenotype = temp[17] ; DisorderCategory = temp[18]
        Disease_category[PrimaryPhenotype] = DisorderCategory
        Disease_list.append(PrimaryPhenotype)
        Category_list.append(DisorderCategory)
    Disease_list = list(set(Disease_list))
    Category_list = list(set(Category_list))

for disease in Disease_list:
    outputdir1 = r'E:\OneDrive\学术\PPI article\data\mapped data\4. DNM_MIPPI_mapping'+'\\'
    outputdir2 = r'_Complex_PPI_mutation_sites.txt'
    disease_name = disease.replace('/','')
    outputdir = outputdir1 +Disease_category[disease]+'\\'+disease_name+outputdir2
    output = open(outputdir,'w+')

    opendir1 = r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping'+'\\'
    opendir2 = r'_DNM_Protein_mutation_mapping.txt'
    disease_name = disease.replace('/','')
    opendir = opendir1+Disease_category[disease]+'\\'+ disease_name + opendir2

    with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\InterfacesHQ_with_detail.txt') as fp:
        output.write('P1'+'\t' + 'P2' + '\t'+ 'Source' +'\t'+ 'P1_symbol'+'\t' +'P2_symbol'+'\t'+'P1_entrez_id'+'\t'+'P2_entrez_id'+'\t'+'P1_mutation_num'+'\t'+'P1_PPI_mutation_num'+'\t'
        +'P2_mutation_num'+'\t'+'P2_PPI_mutation_num'+'\t' + 'P1_PPI_mutations' +'\t'+ 'P2_PPI_mutations'+'\t'+'P1_PPI_mutations_SIFT'+'\t'+'P2_PPI_mutations_SIFT'+'\t'+'P1_PPI_mutations_PolyPhen'+'\t'+'P2_PPI_mutations_PolyPhen'+'\t'
        +'p value' +'\t'+ "P1_pred_disrupting" + '\t' + "P1_pred_decreasing"  + '\t' + "P1_pred_increasing"+ '\t' + "P1_pred_no_effect" + '\t' + "P2_pred_disrupting" + '\t' + "P2_pred_decreasing"
        + '\t' + "P2_pred_increasing"+ '\t' + "P2_pred_no_effect" + '\n')
        text = fp.readline()
        for line in fp:
            P1_mutation_num = 0
            P2_mutation_num = 0
            P1_PPI_mutation_num = 0
            P2_PPI_mutation_num = 0
            interfaces_P1 = []
            interfaces_P2 = []
            P1_PPI_mutations = []
            P2_PPI_mutations = []
            P1_PPI_mutations_SIFT = []
            P2_PPI_mutations_SIFT =[]
            P1_PPI_mutations_PolyPhen = []
            P2_PPI_mutations_PolyPhen =[]
            P1_AAchange_Mutation_on_interface = []
            P2_AAchange_Mutation_on_interface = []
            P1_pred_decreasing = 0
            P1_pred_disrupting = 0
            P1_pred_no_effect = 0
            P1_pred_increasing = 0
            P2_pred_decreasing = 0
            P2_pred_disrupting = 0
            P2_pred_no_effect = 0
            P2_pred_increasing = 0
            pvalue = 1
            temp = line.split('\t')
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

            with open(opendir) as md:
                for row in md:
                    P1_mutation_list = []
                    P2_mutation_list = []
                    P1_SIFT_score = ''
                    P2_SIFT_score = ''
                    P1_Poly_score = ''
                    P2_Poly_score = ''
                    P1_PPI_mutation_num_judge = 0
                    P2_PPI_mutation_num_judge = 0
                    content = row.strip()
                    mapping = content.split('\t')
                    UniprotID = mapping[0]; GeneID = mapping[1] ;Gene_refGene = mapping[2]; PrimaryPhenotype = mapping[4] ; DisorderCategory = mapping[5]; AAchange = mapping[6]
                    SIFT_score = mapping[7]; Polyphen2_HDIV_score = mapping[8] ;CADD_phred = mapping[9]; Protein_Change = mapping[10]

                    if P1 == UniprotID:
                        round = 0
                        P1_symbol = Gene_refGene
                        P1_mutation_num += 1
                        P1_Protein_Change = Protein_Change.strip()
                        P1_Protein_Change_list = P1_Protein_Change.split(',')
                        P1_Protein_Change_list = list(set(P1_Protein_Change_list))
                        for i in P1_Protein_Change_list:
                            P1_mutation_site = re.findall(r'\d+',i)

                            if int(P1_mutation_site[0]) in interfaces_P1:
                                P1_PPI_mutation_num_judge = 1
                                P1_PPI_mutations.append(i)

                        if P1_PPI_mutation_num_judge == 1:
                            P1_AAchange_Mutation_on_interface.append(AAchange)
                            P1_PPI_mutation_num  += 1
                            P1_SIFT_score = SIFT_score
                            P1_Poly_score = Polyphen2_HDIV_score

                            if P1_SIFT_score != '.':
                                if float(P1_SIFT_score) <=0.05:
                                    type_SIFT = 'deleterious'
                                elif float(P1_SIFT_score) >0.5:
                                    type_SIFT = 'tolerated'
                            else:
                                type_SIFT = 'unknown'
                            if P1_Poly_score !='.':
                                if float(P1_Poly_score) >=0.957:
                                    type_Poly = 'probably_damaging'
                                elif  0.453<=float(P1_Poly_score)< 0.957 :
                                    type_Poly = 'possibly_damaging'
                                elif float(P1_Poly_score) < 0.453:
                                    type_Poly = 'benign'
                            else:
                                type_Poly = 'unknown'

                            P1_PPI_mutations_SIFT.append(type_SIFT+'('+str(P1_SIFT_score)+')')
                            P1_PPI_mutations_PolyPhen.append(type_Poly+'('+str(P1_Poly_score)+')')

                    if P2 == UniprotID:
                        round = 0
                        P2_symbol = Gene_refGene
                        P2_mutation_num += 1
                        P2_Protein_Change = Protein_Change.strip()
                        P2_Protein_Change_list = P2_Protein_Change.split(',')
                        P2_Protein_Change_list = list(set(P2_Protein_Change_list))
                        for i in P2_Protein_Change_list:
                            P2_mutation_site = re.findall(r'\d+',i)

                            if int(P2_mutation_site[0]) in interfaces_P2:
                                P2_PPI_mutation_num_judge = 1
                                P2_PPI_mutations.append(i)

                        if P2_PPI_mutation_num_judge == 1:
                            P2_AAchange_Mutation_on_interface.append(AAchange)
                            P2_PPI_mutation_num  += 1
                            P2_SIFT_score = SIFT_score
                            P2_Poly_score = Polyphen2_HDIV_score

                            if P2_SIFT_score != '.':
                                if float(P2_SIFT_score) <=0.05:
                                    type_SIFT = 'deleterious'
                                elif float(P2_SIFT_score) >0.5:
                                    type_SIFT = 'tolerated'
                            else:
                                type_SIFT = 'unknown'
                            if P2_Poly_score !='.':
                                if float(P2_Poly_score) >=0.957:
                                    type_Poly = 'probably_damaging'
                                elif  0.453<=float(P2_Poly_score)< 0.957 :
                                    type_Poly = 'possibly_damaging'
                                elif float(P2_Poly_score) < 0.453:
                                    type_Poly = 'benign'
                            else:
                                type_Poly = 'unknown'

                            P2_PPI_mutations_SIFT.append(type_SIFT+'('+str(P2_SIFT_score)+')')
                            P2_PPI_mutations_PolyPhen.append(type_Poly+'('+str(P2_Poly_score)+')')


                if P1_PPI_mutation_num or P2_PPI_mutation_num != 0:
                    P1_PPI_mutation_num = int(P1_PPI_mutation_num)
                    P2_PPI_mutation_num = int(P2_PPI_mutation_num)
                    P1_mutation_num = int(P1_mutation_num)
                    P2_mutation_num = int(P2_mutation_num)
                    binom_p1 = int(P1_binding_length)/int(P1_protein_length)
                    binom_p2 = int(P2_binding_length)/int(P2_protein_length)

                    p1_value = 1
                    if P1_PPI_mutation_num and P1_mutation_num >0:
                        p1_value = stats.binom_test(P1_PPI_mutation_num, P1_mutation_num, binom_p1 , alternative='greater')
                    p2_value = 1
                    if P2_PPI_mutation_num and P2_mutation_num >0:
                        p2_value = stats.binom_test(P2_PPI_mutation_num, P2_mutation_num, binom_p2 , alternative='greater')
                    if binom_p1 == 0:
                        p1_value = 1
                    if binom_p2 == 0:
                        p2_value =1
                    pvalue = p1_value * p2_value

                    for change in P1_AAchange_Mutation_on_interface:
                        with open(r'E:\OneDrive\学术\PPI article\data\MIPPI\Uniprot_mappinged_MIPPI_on_interface_insiderfiltered.txt') as sd:
                            sd.readline()
                            for text in sd:
                                info = text.strip()
                                detail = info.split('\t')
                                AAchange_ref = detail[0]; SourceUni = detail[1]; TargetUni = detail[3] ; pred_class = detail[5]
                                if (P1 == SourceUni and P2 ==TargetUni) or (P2 == SourceUni and P1 ==TargetUni):
                                    if AAchange_ref in change:
                                        if pred_class == "decreasing":
                                            P1_pred_decreasing += 1
                                        elif pred_class == "disrupting":
                                            P1_pred_disrupting += 1
                                        elif pred_class == "no effect":
                                            P1_pred_no_effect += 1
                                        elif pred_class == "increasing":
                                            P1_pred_increasing += 1
                                        break
                    for change in P2_AAchange_Mutation_on_interface:
                        with open(r'E:\OneDrive\学术\PPI article\data\MIPPI\Uniprot_mappinged_MIPPI_on_interface_insiderfiltered.txt') as sp:
                            sp.readline()
                            for text in sp:
                                info = text.strip()
                                detail = info.split('\t')
                                AAchange_ref = detail[0]; SourceUni = detail[1]; TargetUni = detail[3] ; pred_class = detail[5]
                                if (P2 == SourceUni and P1 ==TargetUni) or (P1 == SourceUni and P2 ==TargetUni):
                                    if AAchange_ref in change:
                                        if pred_class == "decreasing":
                                            P2_pred_decreasing += 1
                                        elif pred_class == "disrupting":
                                            P2_pred_disrupting += 1
                                        elif pred_class == "no effect":
                                            P2_pred_no_effect += 1
                                        elif pred_class == "increasing":
                                            P2_pred_increasing += 1
                                        break

                    output.write(P1 + '\t' + P2 + '\t' + Source + '\t' + P1_symbol + '\t' + P2_symbol + '\t' + P1_entreid + '\t' + P2_entreid  + '\t' +
                    str(P1_mutation_num) + '\t' + str(P1_PPI_mutation_num) + '\t' + str(P2_mutation_num) + '\t' + str(P2_PPI_mutation_num) + '\t' +
                    str(P1_PPI_mutations) + '\t' + str(P2_PPI_mutations) + '\t' + str(P1_PPI_mutations_SIFT) + '\t' + str(P2_PPI_mutations_SIFT) + '\t' +
                    str(P1_PPI_mutations_PolyPhen)  + '\t'+ str(P2_PPI_mutations_PolyPhen) +'\t' +str(pvalue) + '\t' + str(P1_pred_disrupting) + '\t' + str(P1_pred_decreasing)
                    + '\t' + str(P1_pred_increasing)+ '\t' + str(P1_pred_no_effect) + '\t' + str(P2_pred_disrupting) + '\t' + str(P2_pred_decreasing)
                    + '\t' + str(P2_pred_increasing)+ '\t' + str(P2_pred_no_effect) + '\n')

output.close()
print('MIPPI_mutation_sites_mapping by disease is done')
