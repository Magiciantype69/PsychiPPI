import re
import sys
import os
import math

IRES_list = []

output = open(r'E:\OneDrive\学术\PPI article\data\Biogrid\result\IDsurfaceID.txt','w+')

with open('E:\OneDrive\学术\PPI article\data\Biogrid\IDsurfaceID.txt') as root:
    root.readline()
    for line in root:
        ID=line.strip()
        IRES_list.append(ID)
    IRES_list = list(set(IRES_list))

    for id in IRES_list:
        id_PPI = []
        with open(r'E:\学术\PPI project\BIOGRID-ALL\BIOGRID-ALL-High.txt') as wood:
            wood.readline()
            for row in wood:
                row = row.strip()
                row = row.split('\t')
                InteractorA = row[0]; InteractorB = row[1]

                if id == InteractorA:
                    id_PPI.append(InteractorB)
                if id == InteractorB:
                    id_PPI.append(InteractorA)
            id_PPI = list(set(id_PPI))
            output.write(str(id)+'\t'+str(len(id_PPI))+'\n')
output.close()
