import re
import sys
import os
import math

print('DNM_on_interface&MIPPI')
output = open(r'E:\OneDrive\学术\PPI article\data\MIPPI\MIPPI&DNM_only_on_interface.txt','w+')

with open(r'E:\OneDrive\学术\PPI article\data\MIPPI\DNM_on_ppi_interface_filter.txt') as fp:
    fp.readline()
    for line in fp:
        temp = line.split('\t')
        AAchange = temp[6]
        with open(r'E:\OneDrive\学术\PPI article\data\MIPPI\merged_filtered.txt') as md:
            for row in md:
                mippi = row.split('\t')
                mippi_AA = mippi[0]
                mippi_partner = mippi[2]
                mippi_pred = mippi[3]
                mippi_pred_score = mippi[4]
                if mippi_AA in AAchange:
                    output.write(row)
    output.close()
    print('done')
