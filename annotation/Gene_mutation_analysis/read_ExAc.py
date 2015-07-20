
#!/usr/bin/python
import os, glob,sys, subprocess,re,math,random 

## read VCF and use tabix to get mutations in a genomic interval chr:firstpos-lastpos corresponding to a gene 
## update allele frequency info in 'mutations' table/index 

## this can read any VCF file: ExAC or NHLBI exomes, it should be indexed using tabix and bgzip 

## pass mutations hashtable to it and it is updated
def read_vcf(vcffile,chrom,firstpos,lastpos,mutations):

	#################################################################################################
	## run tabix to get mutation data from EXAC VCF file
	#tabix -p vcf /media/drive2/Variant-call-datasets/65Kexomes-DGM.v0.1/ExAC.r0.1.sites.vep.vcf.gz  20:47989510-48099180 > KCNB1_exac.tabix

	tabix_command = ["/home/vbansal/Public/tools/htslib-master/tabix"];
	if vcffile == 'None': tabix_command.append("/media/drive2/Variant-call-datasets/65Kexomes-DGM.v0.1/ExAC.r0.1.sites.vep.vcf.gz");
	else: tabix_command.append(vcffile); 

	tabix_command.append(chrom + ':' + `firstpos` + '-' + `lastpos`); 
	#print tabix_command;
	outfilename = 'EXAC.data.' + chrom + ':' + `firstpos` + '-' + `lastpos`;
	f = open(outfilename,'wb');
	subprocess.call(tabix_command,stdout=f);
	f.close();
	variant_list = []; 

	# 20      47989752        .       G       A,C     17618.44        PASS    AC=3,1;AC_AFR=0,1;AC_AMR=0,0;AC_Adj=3,1;
	File= open(outfilename,'r');
	for line in File:
		var = line.strip().split('\t');
		chrom = var[0]; position = var[1]; ref = var[3]; alleles = var[4].split(',');
		info = var[7].split(';'); 
		if info[0].split('=')[0] == 'AC': AC = info[0].split('=')[1].split(','); ac = 1; 
		else: ac = 0;
		for i in xrange(len(alleles)):
			alt = alleles[i]
			#print chrom,position,ref,alt,int(AC[i]);
			try: 
				m = mutations[(var[0],int(var[1]),var[3],alt)]; 
				m[8]= 1; 
			except KeyError: 
				pass; 
				#print 'notfound';
			if ac ==1: variant_list.append([var[0],int(var[1]),var[3],alt,int(AC[i])]); 
			else: variant_list.append([var[0],int(var[1]),var[3],alt,0]);
		#print var[0:5];
	
	return variant_list; 


