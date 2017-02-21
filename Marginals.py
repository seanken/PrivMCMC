from DataIterator import *;
import numpy as np;
import scipy as sp;


##
##Caculates Marginals
##
class Marginals:
	def __init__(self,filename,snpfile,noise,err,numStep=1000,numSamples=1000,params="",analyser="MAF",n0=-1,n1=-1):
		self.iter=DatasetIterator(filename,snpfile,noise,err,params,analyser,True,n0,n1);
		self.numSteps=numStep;
		self.numSamples=numSamples;
		self.n=self.iter.getN();


	def getMarginals(self,randStart=True,burnIn=100000):
		count=np.asarray([0 for i in range(0,self.n)]);
		print "Burn In!"
		self.iter.burnIn(burnIn);
		for i in range(0,self.numSamples):
			print "Sample "+str(i);
			#self.iter.reInitialize(randStart)
			sampY=np.asarray(self.iter.sample(self.numSteps))
			count=count+sampY;
		count=count/float(self.numSamples);
		return count;

	##
	##One test of convergence: compares high marginals between random and not random start
	##If similar probably converged
	##return graph of numSteps versus difference, or similar
	##
	def checkConverges(self):
		print "To be implemented!"
		return True;

	def setExact(self,exact):
		self.iter.setExact(exact);


