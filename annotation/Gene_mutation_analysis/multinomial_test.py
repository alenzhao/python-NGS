import sys,os, math
from optparse import OptionParser
from numpy import random
from numpy import array

## http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.multinomial.html 
## march 23 2015 code to calculate multinomial p-value
## python multinomial_test.py MLL2.stopmutations /media/drive2/Kabuki_mutations/MUTATION_LIST

## what if we calculate joint probability across multiple studies treating them as independent replicates from binomial distribution.... rather than summing up counts 

## mann whitney test: calculate expected rank sum of 'X' mutations sampled from space of 1852 mutations where each mutation has different probability... 
## 00000100000001001000000000000100010100000000 134 '1' selected from 1852 mutations... we have to calculate whether there is any bias 
## expected rank of mutation = (p1.1 + p2.2 + p3.3 + ...) weighted sum of bernoulli random variables
## http://mathforum.org/kb/message.jspa?messageID=7063154 useful for this special case...
## mean and variance easy to calculate... use chenoff bound for 'X' independent i.i.d. random variables...


################################ code for testing each individual exon for excess of mutation ################

def ncr(n,r):
	ll = 0;
	for i in xrange(min(r,n-r)): ll += math.log(float(n-i)/(i+1),10);
	return ll;

def binomial_pvalue(R,A,e): ## R is total number of mutations, A is observed count and 'e' is mutation rate 
	e1 = math.log(e,10); e2 = math.log(1.0-e,10);
	ll = ncr(R,A) + A*e1 + (R-A)*e2; # pdf of (R,A,e) 
	#print '\n',ll,R,A,e;
	pvlog = ll; pvsumr = ll;  
	for r in xrange(A+1,R+1,1):		
		pvlog += math.log(R-r+1,10) - math.log(r,10) + e1-e2;
		if pvsumr > pvlog: pvsumr += math.log(1.0 + math.pow(10,pvlog-pvsumr),10); 
		else: pvsumr = pvlog + math.log(1.0 + math.pow(10,pvsumr-pvlog),10);

	if A == 0: pvsumr = 0; pvsuml = ll; 
	elif A < R: 
		tempv = pvsumr + math.log(1.0-math.pow(10,ll-pvsumr),10); 
		pvsuml = math.log(1.0-math.pow(10,tempv),10);
	else: pvsuml = 0; 

	return [pvsuml,pvsumr];

	"""
	ll = R*e2; pvlog = ll; 
	for r in xrange(1,A,1):
		pvlog += math.log(R-r+1,10) - math.log(r,10) + e1-e2;
		if pvsum > pvlog: pvsum += math.log(1.0 + math.pow(10,pvlog-pvsum),10); 
		else: pvsum = pvlog + math.log(1.0 + math.pow(10,pvsum-pvlog),10);
	"""


def single_exon_test(pvals,counts):
	bins = len(counts);
	mutations = 0;
	for i in xrange(bins): mutations += counts[i]; 

	for i in xrange(bins):
		p = pvals[i]; muts = counts[i]; ## compare to 1.0 and mutations... 
		pvalue = binomial_pvalue(mutations,muts,p); 
		if pvalue[0]< -2 or pvalue[1] < -2: print 'exon',i+1,p,muts,mutations,pvalue[0],pvalue[1],int(p*mutations);

##########################################################################################################

## calculate multinomial probability for counts and given probability vector 
def multinomial_pdf(counts,probs,n):
	## calculate M(n,a1,a2...ak).p1^a1 p2^a2  pk^ak where a1+a2 + ... + ak = n 
	#print n;
	prob = 0.0; c = float(n); 
	for i in xrange(len(counts)):
		#prob += counts[i]*math.log(probs[i]);
		for j in xrange(counts[i]):
			prob += math.log(probs[i]*c)-math.log(j+1);
			c -= 1; 
	return prob; 

def multinomial_test(pvals,counts):
	bins = len(counts);
	mutations = 0;
	for i in xrange(bins): mutations += counts[i]; 

	C =  random.multinomial(mutations,pvals, size=1); 
	for i in xrange(len(counts)): C[0][i] = counts[i];

	pf = multinomial_pdf(counts,pvals,mutations);
	print C[0],pf,mutations,'REAL_DATA'; 


	pvalue = 0.0; pvalue_exon = 0.0; pvalue_exon10_11 = 0.0;
	for iters in xrange(100000): 
		R = random.multinomial(mutations,pvals, size=1)[0];
		pf_r = multinomial_pdf(R,pvals,mutations);
		if pf_r < pf: pvalue +=1; 
		"""
		print R,pf_r;
		if R[12] > counts[12]: pvalue_exon +=1.0; #print '51+';
		if R[1] + R[2]  <  counts[1] + counts[2]: pvalue_exon10_11 +=1.0; #print '10_11+';
		"""

	print >>sys.stderr, pvalue,iters+1,pvalue/(iters+1),pvalue_exon,pvalue_exon/(iters+1),pvalue_exon10_11,pvalue_exon10_11/(iters+1);



