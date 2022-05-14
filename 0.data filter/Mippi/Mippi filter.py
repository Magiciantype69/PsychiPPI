import re
import sys
import os
import math
import csv

Genename_to_ID ={}

with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\mippi_geneID.txt') as wood:
    for line in wood:
        line = line.strip()
        line = line.split('\t')
        Genename = line[0]
        GeneID = line[1]
        Genename_to_ID[Genename] = GeneID


count = 0
count2 = 0
output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\mippi_filter.txt','w+')
with open(r'E:\学术\PPI project\merged_all.csv') as fp:
    fp_csv = csv.reader(fp)
    next(fp_csv)
    output.write('EntrezId'+'\t'+'Gene symbol'+'\t'+'AAchange'+'\t'+'PartnerID'+'\t'+'Partnersymbol'+'\t'+'prediction'+'\t'+'predict_score'+'\n')
    for row in fp_csv:
        if row[0] == '':
            row[0] = '0'
        if row[12] in Genename_to_ID:
            output.write(str(int(float(row[0])))+'\t'+row[1]+'\t'+row[6]+'\t'+Genename_to_ID[row[12]]+'\t'+row[12]+'\t'+row[20]+'\t'+row[21]+'\n')
            count += 1
        count2 += 1

output.close()
print(count)
print(count2)

count3 =0
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\mippi_filter.txt') as very:
    for line in very:
        count3 += 1
        if count3 > 1234300:
            print(line)
print(count3)
