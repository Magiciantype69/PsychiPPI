import re
import sys
import os
import math

print('DNM_on_interface&MIPPI')
output = open(r'E:\OneDrive\学术\PPI article\data\MIPPI\DNM_only_on_interface.txt','w+')


with open(r'E:\OneDrive\学术\PPI article\data\mapped data\1. DNM_on_interface_mapping\All_DNM_Protein_mutation_mapping.txt') as fp:
    output.write(fp.readline().strip('\n')+'\t'+'Disrupting Partner'+'\t'+'Disrupting count'+'\t'+'Decreasing Partner'+'\t'+'Decreasing count'+'\t'+'Increasing Partner'+'\t'+'Increasing count'+'\t'+'no effect Partner'+'\t'+'no effect count'+'\t'+'mippi count'+'\n')
    for line in fp:
        Disrupting_Partner=[]
        Disrupting_count = 0
        Decreasing_Partner=[]
        Decreasing_count = 0
        Increasing_Partner=[]
        Increasing_count = 0
        no_effect_Partner=[]
        no_effect_count = 0
        temp = line.split('\t')
        AAchange = temp[6]
        mippi_count = 0
        with open(r'E:\OneDrive\学术\PPI article\data\MIPPI\merged_filtered.txt') as md:
            md.readline()
            for row in md:
                mippi = row.split('\t')
                mippi_AA = mippi[0]
                mippi_partner = mippi[2]
                mippi_pred = mippi[3]
                mippi_pred_score = mippi[4]
                if mippi_AA in AAchange:
                    if mippi_pred == "disrupting":
                        Disrupting_count+=1
                        Disrupting_Partner.append(mippi_partner+'('+mippi_pred_score+')')
                        mippi_count +=1
                    elif mippi_pred == "decreasing":
                        Decreasing_count+=1
                        Decreasing_Partner.append(mippi_partner+'('+mippi_pred_score+')')
                        mippi_count +=1
                    elif mippi_pred == "increasing":
                        Increasing_count+=1
                        Increasing_Partner.append(mippi_partner+'('+mippi_pred_score+')')
                        mippi_count +=1
                    elif mippi_pred == "no effect":
                        no_effect_count+=1
                        no_effect_Partner.append(mippi_partner+'('+mippi_pred_score+')')
                        mippi_count +=1
            output.write(line.strip('\n')+'\t'+(','.join(Disrupting_Partner))+'\t'+str(Disrupting_count)+'\t'+(','.join(Decreasing_Partner))+'\t'+str(Decreasing_count)+'\t'+(','.join(Increasing_Partner))+'\t'+str(Increasing_count)+'\t'+(','.join(no_effect_Partner))+'\t'+str(no_effect_count)+'\t'+str(mippi_count)+'\n')
output.write('done')
output.close()
print('done')
