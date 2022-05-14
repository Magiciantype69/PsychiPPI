import re
import sys
import os
import math

output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\Uniprot Protein domain.txt','w+')
with open(r'E:\OneDrive\学术\PPI article\data\other\protein domain.txt') as fp:
    output.write(fp.readline())
    for line in fp:
        domaincontent = []
        domain = []
        text = line.strip('\n')
        text = text.split('\t')
        UniprotID = text[0]
        domaininfo = text[1].split(';')
        for info in domaininfo:
            if 'DOMAIN' in info:
                domaincontent.append(info)
        for info in domaincontent:
            #i = info.strip(r'DOMAIN')
            i = re.findall('\d+',info)
            if i != []:
                for x in range(int(i[0]),int(i[1])+1):
                    domain.append(str(x))
        output.write(UniprotID + '\t' +','.join(domain)+'\n')
output.close()

output2 = open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\Uniprot_Genename_All_domain_sites.txt','w+')
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\Uniprot_Genename_All_PPI_sites.txt') as fd:
    output2.write('UniprotID'+'\t'+'Genename'+'\t'+'Domain_sites'+'\n')
    fd.readline()
    for line in fd:
        PPI_sites = []
        line = line.strip('\n')
        text = line.split('\t')
        if len(text) < 3:
            UniProt = text[0] ; Genename = text[1] ; PPI_site = []
        else:
            UniProt = text[0] ; Genename = text[1] ; PPI_site = text[2].split(',')
        with open(r'E:\OneDrive\学术\PPI article\data\mapped data\0.data filter\Uniprot Protein domain.txt') as fx:
            fx.readline()
            for row in fx:
                content = row.strip('\n')
                content = content.split('\t')
                uniprotdomain = content[0]
                if uniprotdomain ==UniProt:
                    if len(content)<2:
                        domain_list = []
                    else :
                        domain_list = content[1].split(',')
                        for i in PPI_site:
                            if i in domain_list:
                                domain_list.remove(i)
                    break
        output2.write(UniProt+'\t'+Genename +'\t'+','.join(domain_list)+'\n')
output2.close()
