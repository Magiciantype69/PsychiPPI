import re
import sys
import os
import math
from Bio import SwissProt
import urllib.parse
import urllib.request
from pyfaidx import Fasta

#mapping Genename with UniprotID

#Matching UniprotID with it's protein lengh from uniprot.fasta
def extract_id(header):
    return header.split('|')[1]
sequences = Fasta(r'C:\Users\Administrator\Desktop\uniprot_sprot.fasta\uniprot_sprot.fasta', key_function=extract_id)

UniprotID_match_Genename = {}
UniprotID_match_protein_length = {}
UniprotID_match_entreid = {}
with open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\other\All_UNItoGENE.txt') as fp:
    fp.readline()
    for line in fp:
        line = line.strip()
        temp = line.split("\t")
        Genename = temp[1]
        Uniprotid = temp[0]
        try:
            seq = len(sequences[Uniprotid])
            UniprotID_match_protein_length[Uniprotid] = seq
        except:
            pass

        UniprotID_match_Genename[Uniprotid] = Genename

#Matching the rest of error Uniprot from Uniprot website
with open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\other\error.txt') as fd:
    for line in fd:
        line = line.strip()
        temp = line.split("\t")
        Uniprotid = temp[0]
        error_length = temp[1]
        UniprotID_match_protein_length[Uniprotid] = error_length

with open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\other\entreID.txt') as fp:
    fp.readline()
    for line in fp:
        line = line.strip()
        temp= line.split('\t')
        Uniprotid = temp[0]
        entreid = temp[1]
        UniprotID_match_entreid[Uniprotid] =entreid

#sort each ppi interface
InterfacesHQ_with_detail = open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\mapped data\part1_Uniprot_mapping\InterfacesHQ_with_detail.txt','w+')
InterfacesHQ_with_detail.write('P1'+'\t'+'P2'+'\t'+'Source'+'\t'+'P1_symbol'+'\t'+'P2_symbol'+'\t'+'P1_entreid'+'\t'+'P2_entreid'+'\t'+'P1_PPI_residue'+'\t'+'P2_PPI_residue'+'\t'+'P1_protein_length'+'\t'+'P2_protein_length'+'\t'+'P1_binding_length'+'\t'+'P2_bind_length'+'\n')
with open(r'C:\Users\Administrator\OneDrive\学术\PPI project\data\PPI interface data\H_sapiens_interfacesHQ.txt') as fp:
    fp.readline()
    for line in fp:
        line = line.strip()
        temp = line.split('\t')
        Source = temp[2]
        P1 = temp[0]
        P1_IRES = temp[3]
        P1_PPI_site = []
        P2 = temp[1]
        P2_IRES = temp[4]
        P2_PPI_site = []

        try:
            if P1 or P2 in UniprotID_match_Genename:
                if P1_IRES!='[]':
                    P1_IRES = P1_IRES[1:-1].split(',')
                    for IRES in P1_IRES:
                        if "-" in IRES:
                            IRES = IRES.split('-')
                            for i in range(int(IRES[0]),int(IRES[1])+1):
                                P1_PPI_site.append(i)
                        else:
                            P1_PPI_site.append(int(IRES))
                if P2_IRES!='[]':
                    P2_IRES = P2_IRES[1:-1].split(',')
                    for IRES in P2_IRES:
                        if "-" in IRES:
                            IRES = IRES.split('-')
                            for i in range(int(IRES[0]),int(IRES[1])+1):
                                P2_PPI_site.append(i)
                        else:
                            P2_PPI_site.append(int(IRES))
            P1_bind = str(len(P1_PPI_site))
            P2_bind = str(len(P2_PPI_site))
            P1_length = str(UniprotID_match_protein_length[P1])
            P2_length = str(UniprotID_match_protein_length[P2])
            P1_symbol = UniprotID_match_Genename[P1]
            P2_symbol = UniprotID_match_Genename[P2]
            P1_entreid = UniprotID_match_entreid[P1]
            P2_entreid = UniprotID_match_entreid[P2]
            if len(P1_PPI_site) or len(P2_PPI_site) !=0:
                InterfacesHQ_with_detail.write(P1+'\t'+P2+'\t'+Source+'\t'+
                P1_symbol+'\t'+P2_symbol+'\t'+P1_entreid+'\t'+P2_entreid+'\t'+str(P1_PPI_site)+'\t'+str(P2_PPI_site)+'\t'+
                P1_length+'\t'+P2_length+'\t'+P1_bind+'\t'+P2_bind+'\n')
        except:
            pass

InterfacesHQ_with_detail.close()
