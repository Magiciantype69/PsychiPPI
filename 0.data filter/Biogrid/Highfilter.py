import re
import sys
import os
import math



output = open(r'E:\学术\PPI project\BIOGRID-ALL\BIOGRID-ALL-High.txt','w+')

count = 0
count2 =0
with open(r'E:\学术\PPI project\BIOGRID-ALL\BIOGRID-ALL-4.4.205.tab2.txt') as wood:
    title = wood.readline()
    title = title.strip()
    line = title.split('\t')
    output.write(line[1]+'\t'+line[2]+'\t'+line[7]+'\t'+line[8]+'\t'+line[17]+'\n')
    for line in wood:
        line = line.strip()
        line = line.split('\t')
        InteractorA = line[1];InteractorB = line[2];Asymbol = line[7];Bsymbol = line[8]; Throughput = line[17]; Organism =line[16]
        count2 +=1
        if Organism == '9606':
            if 'High' in Throughput:
                output.write(line[1]+'\t'+line[2]+'\t'+line[7]+'\t'+line[8]+'\t'+line[17]+'\n')
                count += 1

output.close()
print(count)
print(count2)
