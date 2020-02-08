import numpy as np
import pandas as pd

# making list of chromosomes in vcf

with open('1000G_phase3_v4_20130502.sites.vcf') as snp:
    c = set()
    for lin in snp:
        lin = lin.rstrip().split()
        if '#' not in lin[0]:
            c.add(lin[0])
chrom = []
for elem in c:
    chrom.append(elem)
chrom = [int(i) for i in chrom]
chrom.sort()

# function that makes dataframe in csv consisting of next columns:
# chromosome, start and end of interval, mean, sd, and length of each interval for every 1 billion of
# positions in one chromosome


def stat(chromosome):
    d = dict()
    d['chrom'] = []
    d['pos'] = []
    d['af'] = []
    with open('1000G_phase3_v4_20130502.sites.vcf') as f:
        for line in f:
            line = line.rstrip().split()
            if '#' not in line[0] and int(line[0]) == chromosome:
                d['chrom'].append(line[0])
                d['pos'].append(int(line[1]))
                AF = line[7].split(';')[1].split('=')[1]
                AF = np.mean([float(i) for i in AF.split(',')])
                d['af'].append(AF)
    df = pd.DataFrame(d)
    n = 1
    snp = dict()
    snp['chromosome'] = []
    snp['start'] = []
    snp['end'] = []
    snp['mean'] = []
    snp['sd'] = []
    snp['length'] = []
    while n < df['pos'].iloc[-1]:
        XAF = df[(n <= df.pos) & (df.pos < n + 999999)]['af']
        if XAF.shape[0] == 0:
            snp['chromosome'].append(df.chrom[0])
            snp['start'].append(n)
            snp['end'].append(n + 999999)
            snp['mean'].append(0)
            snp['sd'].append(0)
            snp['length'].append(0)
            n += 1000000
        else:
            AF_mean = XAF.mean()
            AF_sd = XAF.std()
            AF_len = XAF.shape[0]
            snp['chromosome'].append(df.chrom[0])
            snp['start'].append(n)
            snp['end'].append(n + 999999)
            snp['mean'].append(AF_mean)
            snp['sd'].append(AF_sd)
            snp['length'].append(AF_len)
            # print(f'chr{df.chrom[0]}: ({n} - {n + 999999}) mean = {AF_mean}, sd = {AF_sd}, number = {AF_len}')
            n += 1000000

    data = pd.DataFrame(snp)
    data.to_csv('snp.csv', header=False, sep='\t', encoding='utf-8', index=False, mode='a')
    return

# arranging dataframe for every chromosome in csv


for elem in chrom:
    stat(elem)

