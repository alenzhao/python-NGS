#! /usr/bin/env python

import sys, os, glob, string, subprocess,time, math, re, compiler, random

#/opt/biotools/plink/bin/plink --bfile ASN_17 --make-bed --keep 1198samples --out HapMap3_1198_249075
#python ~/CODE/JOINTCODE-coral/ancestry/runancestry.py --plink European.100sim -f allpops.frq.249075 -o sim -p 2 --path ~/CODE/JOINTCODE-coral/ancestry/
# /opt/biotools/plink/bin/plink --bfile HapMap3_1198_249075 --merge Asian.100sim.ped Asian.100sim.map --make-bed --out HapMap3_Asian_sim
# ../HapMap-AF/HapMap3-admixture-jan2014/admixture_linux-1.23/admixture HapMap3_Asian_sim.bed 11 --supervised

## read allele frequencies from file and generate pl
## argument 2 is YRI:0.2,CEU:0.4,CHB:0.4 
## python simulateadmixed.py allpops.frq.249075 CEU:0.5,TSI:0.5 100 > European.100sim.ped

def read_af(affile):
        aflist = []
        File = open(affile,'r');
        for line in File:
                snp = line.strip().split();
                aflist.append(snp);
        File.close();
        return aflist;

def sample_admix(aflist,admixvec,samples):

        pops = {};
        for i in xrange(5,len(aflist[0])): pops[aflist[0][i]] = i;

        P = [pops[x.split(':')[0]] for x in admixvec.split(',')]; ## population IDS index
        Q = [float(x.split(':')[1]) for x in admixvec.split(',')]; ## population admixture proportions..
        print >>sys.stderr, P,Q,aflist[0];# print Q;

        for s in xrange(samples):
                genotypes = []; LL = 0; LL1 =0;
                print 'SIM_'+admixvec,'SIM_'+`s`,0,0,2,-9,
                for i in range(1,len(aflist)):
                        maf = 0.0;
                        for p in xrange(len(P)): maf += float(aflist[i][P[p]])*Q[p];
                        r = random.random();
                        if r <= maf: allele0 = aflist[i][3];
                        else: allele0 = aflist[i][4];
                        #print aflist[i][5],aflist[i][6],
                        #print r,maf,allele0,

                        r = random.random();
                        if r <= maf: allele1 = aflist[i][3];
                        else: allele1 = aflist[i][4];
                        #print r,maf,allele1;


                        if allele0 == allele1 and allele0 == aflist[i][3]: LL += 2*math.log(maf);
                        elif allele0 == allele1 and allele0 == aflist[i][4]: LL += 2*math.log(1.0-maf);
                        else: LL += math.log(maf) + math.log(1.0-maf) + math.log(2.0);

                        maf = 0.15*float(aflist[i][5]) + 0.85*float(aflist[i][6]);
                        if maf < 0.000001: maf = 0.000001;
                        elif maf > 1.0-0.000001: maf = 1.0-0.000001;
                        if allele0 == allele1 and allele0 == aflist[i][3]: LL1 += 2*math.log(maf);
                        elif allele0 == allele1 and allele0 == aflist[i][4]: LL1 += 2*math.log(1.0-maf);
                        else: LL1 += math.log(maf) + math.log(1.0-maf) + math.log(2.0);

                        print allele0,allele1,
                print;
                print >>sys.stderr, "LL",LL,LL1;


def simulate_AF_matrix(snps,F):
        print >>sys.stdout, "##chrom position rsid A1 A2 CEU TSI";
        for i in xrange(snps):
                p = random.random(); # allele frequency 
                #while p < 0.1 or p > 0.9: p = random.random();
                alpha = (1.0-F)*p; alpha /= F; beta = (1.0-F)*(1.0-p); beta /= F;
                p1 = random.betavariate(alpha,beta);
                p2 = random.betavariate(alpha,beta);
                position = i*100+100;
                print >>sys.stdout, "chr1",position,'rs_' + `position`,'A C',p1,p2;

random.seed();

if len(sys.argv) > 3:
        affile = sys.argv[1]; aflist = read_af(affile);
        sample_admix(aflist,sys.argv[2],int(sys.argv[3]));
else: simulate_AF_matrix(int(sys.argv[1]),float(sys.argv[2]));

