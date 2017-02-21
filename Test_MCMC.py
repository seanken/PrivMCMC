from DataIterator import *;
from pysnptools.snpreader import Bed;
import numpy as np;
import scipy as sp;
from Marginals import *;


def firstTest():
	filename="../GWAS/cleanMAF";
	snpfile="SNP_File.txt";
	numStep=1000;
	numSamples=100;
	n0=1000
	n1=50
	randStart=False;
	numSNP=50;
	X=Bed(filename);
	snps=[X.sid[i] for i in range(0,numSNP)]
	fil=open(snpfile,"w");
	for s in snps:
		fil.write(s+"\n");
	fil.close()
	params="";
	noise=.01;
	err=-1.0;#.05;
	analyser="MAF"
	"""
	iter=DatasetIterator(filename,snpfile,noise,err,params,analyser)
	iter.step()
	print "Took step!"
	T=iter.sample(numSteps=100);
	print T;
	"""
	lst=[];
	for numSNP in [10,20,30,40,50]:
		snps=[X.sid[i] for i in range(0,numSNP)]
		fil=open(snpfile,"w");
		for s in snps:
			fil.write(s+"\n");
		fil.close()
		marg=Marginals(filename,snpfile,noise,err,numStep,numSamples,params,analyser,n0,n1)
		marginals=marg.getMarginals(True);
		#marginals2=marg.getMarginals(False);
		mx=max(marginals)
		#print min(marginals)
		#print max(marginals2)
		#print min(marginals2)
		lst.append(mx)
	print lst;



if __name__=="__main__":
	firstTest();

