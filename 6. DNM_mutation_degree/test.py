import re
import sys
import os
import math

print('DNM_mapping_partner(SCZ) for degree counting')

test = []
with open(r'E:\OneDrive\学术\PPI article\data\mapped data\6. DNM_mapping_for_degree\SCZ.txt') as root:
    title = root.readline()
    for line in root:
        test.append(line)

test = list(set(test))
output = open(r'E:\OneDrive\学术\PPI article\data\mapped data\6. DNM_mapping_for_degree\SCZ_only.txt','w+')

output.write(title)

for i in test:
    output.write(i)

output.close()
