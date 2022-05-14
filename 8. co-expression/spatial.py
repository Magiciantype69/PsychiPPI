import re
import sys
import os
import math
from numpy import *

output=  open(r'E:\OneDrive\学术\博士论文\fig\Part 2\co-expression\SPTBN2 matrix.txt','w+')
p2 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p3 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p4 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p5 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p6 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p7 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p8 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p9 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p10 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p11 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p12 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p13 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p14 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
p15 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

with open(r'E:\OneDrive\学术\博士论文\fig\Part 2\co-expression\SPTBN2.txt') as root:
    root.readline()
    for line in root:
        text =  line.strip()
        text = text.split('\t')
        SI = float(text[0]); region = text[1]; period = int(text[2])

        if period == 2:
            if region == 'OFC':
                p2[0].append(SI)
            elif region == 'DFC':
                p2[1].append(SI)
            elif region == 'VFC':
                p2[2].append(SI)
            elif region == 'MFC':
                p2[3].append(SI)
            elif region == 'M1C':
                p2[4].append(SI)
            elif region == 'S1C':
                p2[5].append(SI)
            elif region == 'IPC':
                p2[6].append(SI)
            elif region == 'A1C':
                p2[7].append(SI)
            elif region == 'STC':
                p2[8].append(SI)
            elif region == 'ITC':
                p2[9].append(SI)
            elif region == 'V1C':
                p2[10].append(SI)
            elif region == 'HIP':
                p2[11].append(SI)
            elif region == 'AMY':
                p2[12].append(SI)
            elif region == 'STR':
                p2[13].append(SI)
            elif region == 'MD':
                p2[14].append(SI)
            elif region == 'CBC':
                p2[15].append(SI)
        if period == 3:
            if region == 'OFC':
                p3[0].append(SI)
            elif region == 'DFC':
                p3[1].append(SI)
            elif region == 'VFC':
                p3[2].append(SI)
            elif region == 'MFC':
                p3[3].append(SI)
            elif region == 'M1C':
                p3[4].append(SI)
            elif region == 'S1C':
                p3[5].append(SI)
            elif region == 'IPC':
                p3[6].append(SI)
            elif region == 'A1C':
                p3[7].append(SI)
            elif region == 'STC':
                p3[8].append(SI)
            elif region == 'ITC':
                p3[9].append(SI)
            elif region == 'V1C':
                p3[10].append(SI)
            elif region == 'HIP':
                p3[11].append(SI)
            elif region == 'AMY':
                p3[12].append(SI)
            elif region == 'STR':
                p3[13].append(SI)
            elif region == 'MD':
                p3[14].append(SI)
            elif region == 'CBC':
                p3[15].append(SI)
        if period == 4:
            if region == 'OFC':
                p4[0].append(SI)
            elif region == 'DFC':
                p4[1].append(SI)
            elif region == 'VFC':
                p4[2].append(SI)
            elif region == 'MFC':
                p4[3].append(SI)
            elif region == 'M1C':
                p4[4].append(SI)
            elif region == 'S1C':
                p4[5].append(SI)
            elif region == 'IPC':
                p4[6].append(SI)
            elif region == 'A1C':
                p4[7].append(SI)
            elif region == 'STC':
                p4[8].append(SI)
            elif region == 'ITC':
                p4[9].append(SI)
            elif region == 'V1C':
                p4[10].append(SI)
            elif region == 'HIP':
                p4[11].append(SI)
            elif region == 'AMY':
                p4[12].append(SI)
            elif region == 'STR':
                p4[13].append(SI)
            elif region == 'MD':
                p4[14].append(SI)
            elif region == 'CBC':
                p4[15].append(SI)
        if period == 5:
            if region == 'OFC':
                p5[0].append(SI)
            elif region == 'DFC':
                p5[1].append(SI)
            elif region == 'VFC':
                p5[2].append(SI)
            elif region == 'MFC':
                p5[3].append(SI)
            elif region == 'M1C':
                p5[4].append(SI)
            elif region == 'S1C':
                p5[5].append(SI)
            elif region == 'IPC':
                p5[6].append(SI)
            elif region == 'A1C':
                p5[7].append(SI)
            elif region == 'STC':
                p5[8].append(SI)
            elif region == 'ITC':
                p5[9].append(SI)
            elif region == 'V1C':
                p5[10].append(SI)
            elif region == 'HIP':
                p5[11].append(SI)
            elif region == 'AMY':
                p5[12].append(SI)
            elif region == 'STR':
                p5[13].append(SI)
            elif region == 'MD':
                p5[14].append(SI)
            elif region == 'CBC':
                p5[15].append(SI)
        if period == 6:
            if region == 'OFC':
                p6[0].append(SI)
            elif region == 'DFC':
                p6[1].append(SI)
            elif region == 'VFC':
                p6[2].append(SI)
            elif region == 'MFC':
                p6[3].append(SI)
            elif region == 'M1C':
                p6[4].append(SI)
            elif region == 'S1C':
                p6[5].append(SI)
            elif region == 'IPC':
                p6[6].append(SI)
            elif region == 'A1C':
                p6[7].append(SI)
            elif region == 'STC':
                p6[8].append(SI)
            elif region == 'ITC':
                p6[9].append(SI)
            elif region == 'V1C':
                p6[10].append(SI)
            elif region == 'HIP':
                p6[11].append(SI)
            elif region == 'AMY':
                p6[12].append(SI)
            elif region == 'STR':
                p6[13].append(SI)
            elif region == 'MD':
                p6[14].append(SI)
            elif region == 'CBC':
                p6[15].append(SI)
        if period == 7:
            if region == 'OFC':
                p7[0].append(SI)
            elif region == 'DFC':
                p7[1].append(SI)
            elif region == 'VFC':
                p7[2].append(SI)
            elif region == 'MFC':
                p7[3].append(SI)
            elif region == 'M1C':
                p7[4].append(SI)
            elif region == 'S1C':
                p7[5].append(SI)
            elif region == 'IPC':
                p7[6].append(SI)
            elif region == 'A1C':
                p7[7].append(SI)
            elif region == 'STC':
                p7[8].append(SI)
            elif region == 'ITC':
                p7[9].append(SI)
            elif region == 'V1C':
                p7[10].append(SI)
            elif region == 'HIP':
                p7[11].append(SI)
            elif region == 'AMY':
                p7[12].append(SI)
            elif region == 'STR':
                p7[13].append(SI)
            elif region == 'MD':
                p7[14].append(SI)
            elif region == 'CBC':
                p7[15].append(SI)
        if period == 8:
            if region == 'OFC':
                p8[0].append(SI)
            elif region == 'DFC':
                p8[1].append(SI)
            elif region == 'VFC':
                p8[2].append(SI)
            elif region == 'MFC':
                p8[3].append(SI)
            elif region == 'M1C':
                p8[4].append(SI)
            elif region == 'S1C':
                p8[5].append(SI)
            elif region == 'IPC':
                p8[6].append(SI)
            elif region == 'A1C':
                p8[7].append(SI)
            elif region == 'STC':
                p8[8].append(SI)
            elif region == 'ITC':
                p8[9].append(SI)
            elif region == 'V1C':
                p8[10].append(SI)
            elif region == 'HIP':
                p8[11].append(SI)
            elif region == 'AMY':
                p8[12].append(SI)
            elif region == 'STR':
                p8[13].append(SI)
            elif region == 'MD':
                p8[14].append(SI)
            elif region == 'CBC':
                p8[15].append(SI)
        if period == 9:
            if region == 'OFC':
                p9[0].append(SI)
            elif region == 'DFC':
                p9[1].append(SI)
            elif region == 'VFC':
                p9[2].append(SI)
            elif region == 'MFC':
                p9[3].append(SI)
            elif region == 'M1C':
                p9[4].append(SI)
            elif region == 'S1C':
                p9[5].append(SI)
            elif region == 'IPC':
                p9[6].append(SI)
            elif region == 'A1C':
                p9[7].append(SI)
            elif region == 'STC':
                p9[8].append(SI)
            elif region == 'ITC':
                p9[9].append(SI)
            elif region == 'V1C':
                p9[10].append(SI)
            elif region == 'HIP':
                p9[11].append(SI)
            elif region == 'AMY':
                p9[12].append(SI)
            elif region == 'STR':
                p9[13].append(SI)
            elif region == 'MD':
                p9[14].append(SI)
            elif region == 'CBC':
                p9[15].append(SI)
        if period == 10:
            if region == 'OFC':
                p10[0].append(SI)
            elif region == 'DFC':
                p10[1].append(SI)
            elif region == 'VFC':
                p10[2].append(SI)
            elif region == 'MFC':
                p10[3].append(SI)
            elif region == 'M1C':
                p10[4].append(SI)
            elif region == 'S1C':
                p10[5].append(SI)
            elif region == 'IPC':
                p10[6].append(SI)
            elif region == 'A1C':
                p10[7].append(SI)
            elif region == 'STC':
                p10[8].append(SI)
            elif region == 'ITC':
                p10[9].append(SI)
            elif region == 'V1C':
                p10[10].append(SI)
            elif region == 'HIP':
                p10[11].append(SI)
            elif region == 'AMY':
                p10[12].append(SI)
            elif region == 'STR':
                p10[13].append(SI)
            elif region == 'MD':
                p10[14].append(SI)
            elif region == 'CBC':
                p10[15].append(SI)
        if period == 11:
            if region == 'OFC':
                p11[0].append(SI)
            elif region == 'DFC':
                p11[1].append(SI)
            elif region == 'VFC':
                p11[2].append(SI)
            elif region == 'MFC':
                p11[3].append(SI)
            elif region == 'M1C':
                p11[4].append(SI)
            elif region == 'S1C':
                p11[5].append(SI)
            elif region == 'IPC':
                p11[6].append(SI)
            elif region == 'A1C':
                p11[7].append(SI)
            elif region == 'STC':
                p11[8].append(SI)
            elif region == 'ITC':
                p11[9].append(SI)
            elif region == 'V1C':
                p11[10].append(SI)
            elif region == 'HIP':
                p11[11].append(SI)
            elif region == 'AMY':
                p11[12].append(SI)
            elif region == 'STR':
                p11[13].append(SI)
            elif region == 'MD':
                p11[14].append(SI)
            elif region == 'CBC':
                p11[15].append(SI)
        if period == 12:
            if region == 'OFC':
                p12[0].append(SI)
            elif region == 'DFC':
                p12[1].append(SI)
            elif region == 'VFC':
                p12[2].append(SI)
            elif region == 'MFC':
                p12[3].append(SI)
            elif region == 'M1C':
                p12[4].append(SI)
            elif region == 'S1C':
                p12[5].append(SI)
            elif region == 'IPC':
                p12[6].append(SI)
            elif region == 'A1C':
                p12[7].append(SI)
            elif region == 'STC':
                p12[8].append(SI)
            elif region == 'ITC':
                p12[9].append(SI)
            elif region == 'V1C':
                p12[10].append(SI)
            elif region == 'HIP':
                p12[11].append(SI)
            elif region == 'AMY':
                p12[12].append(SI)
            elif region == 'STR':
                p12[13].append(SI)
            elif region == 'MD':
                p12[14].append(SI)
            elif region == 'CBC':
                p12[15].append(SI)
        if period == 13:
            if region == 'OFC':
                p13[0].append(SI)
            elif region == 'DFC':
                p13[1].append(SI)
            elif region == 'VFC':
                p13[2].append(SI)
            elif region == 'MFC':
                p13[3].append(SI)
            elif region == 'M1C':
                p13[4].append(SI)
            elif region == 'S1C':
                p13[5].append(SI)
            elif region == 'IPC':
                p13[6].append(SI)
            elif region == 'A1C':
                p13[7].append(SI)
            elif region == 'STC':
                p13[8].append(SI)
            elif region == 'ITC':
                p13[9].append(SI)
            elif region == 'V1C':
                p13[10].append(SI)
            elif region == 'HIP':
                p13[11].append(SI)
            elif region == 'AMY':
                p13[12].append(SI)
            elif region == 'STR':
                p13[13].append(SI)
            elif region == 'MD':
                p13[14].append(SI)
            elif region == 'CBC':
                p13[15].append(SI)
        if period == 14:
            if region == 'OFC':
                p14[0].append(SI)
            elif region == 'DFC':
                p14[1].append(SI)
            elif region == 'VFC':
                p14[2].append(SI)
            elif region == 'MFC':
                p14[3].append(SI)
            elif region == 'M1C':
                p14[4].append(SI)
            elif region == 'S1C':
                p14[5].append(SI)
            elif region == 'IPC':
                p14[6].append(SI)
            elif region == 'A1C':
                p14[7].append(SI)
            elif region == 'STC':
                p14[8].append(SI)
            elif region == 'ITC':
                p14[9].append(SI)
            elif region == 'V1C':
                p14[10].append(SI)
            elif region == 'HIP':
                p14[11].append(SI)
            elif region == 'AMY':
                p14[12].append(SI)
            elif region == 'STR':
                p14[13].append(SI)
            elif region == 'MD':
                p14[14].append(SI)
            elif region == 'CBC':
                p14[15].append(SI)
        if period == 15:
            if region == 'OFC':
                p15[0].append(SI)
            elif region == 'DFC':
                p15[1].append(SI)
            elif region == 'VFC':
                p15[2].append(SI)
            elif region == 'MFC':
                p15[3].append(SI)
            elif region == 'M1C':
                p15[4].append(SI)
            elif region == 'S1C':
                p15[5].append(SI)
            elif region == 'IPC':
                p15[6].append(SI)
            elif region == 'A1C':
                p15[7].append(SI)
            elif region == 'STC':
                p15[8].append(SI)
            elif region == 'ITC':
                p15[9].append(SI)
            elif region == 'V1C':
                p15[10].append(SI)
            elif region == 'HIP':
                p15[11].append(SI)
            elif region == 'AMY':
                p15[12].append(SI)
            elif region == 'STR':
                p15[13].append(SI)
            elif region == 'MD':
                p15[14].append(SI)
            elif region == 'CBC':
                p15[15].append(SI)

for i in p2:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p3:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p4:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p5:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p6:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p7:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p8:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p9:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p10:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p11:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p12:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p13:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p14:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')
for i in p15:
    if len(i) !=0:
        output.write(str(mean(i))+'\t')
    else:
        output.write('nah'+'\t')
output.write('\n')

output.close()
