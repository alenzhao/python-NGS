

def get_mutation_rates():
	## read file 
	PATH="/home/vbansal/CODE/JOINTCODE-coral/FET_contingencytable/SingleIndividual-pvalue/"

	mut_table = {}; bases = ['A','C','G','T']
	File = open(PATH + "trimer_mu_matrix.lis",'r');
	for line in File:
		if line[0] == '#': continue;
		trimer = line.split();
		for i in xrange(4): mut_table[(trimer[0],bases[i])] = float(trimer[i+1]);
	File.close();

	#for trimer in mut_table.keys(): print trimer,mut_table[trimer]
	return mut_table; 



