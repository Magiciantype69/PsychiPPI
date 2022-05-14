import re
import sys
import os
import math
from Bio import SwissProt
import urllib.parse
import urllib.request
from pyfaidx import Fasta

###update the Uniprot ID via uniprot_sprot database
UniProt_ID_update = {}
for record in SwissProt.parse(open(r'C:\Users\Administrator\Desktop\uniprot_sprot.dat\uniprot_sprot.dat')):
    for accession in record.accessions:
        UniProt_ID_update[accession] = record.accessions[0]
print("UniProt ID update is done")

### match each protein's possible interface residue sites with their Uniprot ID
UniProt_interfaces={}
with open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\PPI interface data\H_sapiens_interfacesHQ.txt') as fp:
    fp.readline()
    for line in fp:
        line = line.strip()
        temp = line.split("\t")
        P1 = temp[0]
        P1_IRES = temp[3]
        P2 = temp[1]
        P2_IRES = temp[4]
        if P1 in UniProt_ID_update:
            P1 = UniProt_ID_update[P1]
        if P2 in UniProt_ID_update:
            P2 = UniProt_ID_update[P2]
        if P1 not in UniProt_interfaces:
            UniProt_interfaces[P1] = set()
        if P2 not in UniProt_interfaces:
            UniProt_interfaces[P2] = set()
        if P1_IRES!='[]':
            P1_IRES = P1_IRES[1:-1].split(',')
            for IRES in P1_IRES:
                if "-" in IRES:
                    IRES = IRES.split('-')
                    for i in range(int(IRES[0]),int(IRES[1])+1):
                        UniProt_interfaces[P1].add(str(i))
                else:
                    UniProt_interfaces[P1].add(IRES)
        if P2_IRES!='[]':
            P2_IRES = P2_IRES[1:-1].split(',')
            for IRES in P2_IRES:
                if "-" in IRES:
                    IRES = IRES.split('-')
                    for i in range(int(IRES[0]),int(IRES[1])+1):
                        UniProt_interfaces[P2].add(str(i))
                else:
                    UniProt_interfaces[P2].add(IRES)

PPI_sites = open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\other\Each_protein_all_possible_PPI_sites.txt','w+')
for Protein in UniProt_interfaces:
    PPI_sites.write(Protein + '\t' + ','.join(UniProt_interfaces[Protein]) + '\n')
PPI_sites.close()
print("Uniprot&interface sites matching has been completed successfully")

###Matching Uniprot ID with it's Genename and interface residue
Uniprot_ID_name = ''
Protein_all_PPI_site = {}
with open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\other\Each_protein_all_possible_PPI_sites.txt') as fp:
    for line in fp:
        line = line.strip()
        temp = line.split("\t")
        interface_UniProt_ID= temp[0]
        Uniprot_ID_name = Uniprot_ID_name + interface_UniProt_ID + ' '


url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC+ID',
'to': 'GENENAME',
'format': 'tab',
'query': Uniprot_ID_name
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
  response = f.read()
output = open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\other\UNItoGENE.txt','w+')
output.write(response.decode('utf-8'))
output.close()

#combine the Genename with it's own UniprotID and all possible PPI sites
UniprotID_match_Genename = {}
with open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\other\UNItoGENE.txt') as fp:
    #fp.readline()
    for line in fp:
        line = line.strip()
        temp = line.split("\t")
        Genename1 = temp[1]
        Uniprotid1 = temp[0]
        UniprotID_match_Genename[Uniprotid1] = Genename1

with open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\other\All_UNItoGENE.txt') as fp:
    #fp.readline()
    for line in fp:
        line = line.strip()
        temp = line.split("\t")
        Genename1 = temp[1]
        Uniprotid1 = temp[0]
        UniprotID_match_Genename[Uniprotid1] = Genename1

output = open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\mapped data\part1_Uniprot_mapping\Uniprot_Genename_All_PPI_sites.txt','w+')
output.write('UniprotID'+'\t'+'Genename'+'\t'+'Protein_length'+'\t'+'PPI length'+'\n')
with open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\other\Each_protein_all_possible_PPI_sites.txt') as fp:
    #fp.readline()
    for line in fp:
        line = line.strip()
        temp = line.split("\t")
        Uniprotid = temp[0]
        if Uniprotid not in UniprotID_match_Genename:
            UniprotID_match_Genename[Uniprotid]=''
        if len(temp) == 1:
            output.write((Uniprotid) + '\t' + UniprotID_match_Genename[Uniprotid] +'\n')
        else:
            output.write((Uniprotid) + '\t' +  UniprotID_match_Genename[Uniprotid]+ '\t'+ temp[1]+'\n')
output.close()

print('Genename & UniprotID and corresponding PPI interface sites has been matched')
