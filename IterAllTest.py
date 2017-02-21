import numpy as np;
import scipy as sp;
##from Marginals import *;
import math;
from Analysis import *;



class ItAll:
	def __init__(self,n1,n,filename,snpfile,noise):
		self.y=[0 for i in range(0,n)]
		self.n=n;
		self.n1=n1;
		self.y[:n1]=[1 for i in range(0,n1)]
		self.isDone=False
		self.cases=[i for i in range(0,n1)]

		self.analysis=MAF_Analysis(filename,snpfile);
		self.noisyMAF=self.analysis.getNoisyAnalysis(noise,err=-1)
		self.noise=noise;
		self.m=len(self.noisyMAF);

	def  next(self):
		if self.isDone:
			return;
		I=min([i for i in range(0,self.n1) if self.cases[i]==self.n-1 or self.y[self.cases[i]+1]==0]);
		if self.cases[I]==self.n-1:
			self.isDone=True;
			return;
		self.cases[I]=self.cases[I]+1;
		self.cases[:I]=[i for i in range(0,I)]

	
		self.y=[0 for i in range(0,self.n)]
		for i in self.cases:
			self.y[i]=1;
	
	
	
	def is_Done(self):
		return self.isDone;


	def getY(self):
		return np.array([i for i in self.y])



	def getScore(self):
		maf=self.analysis.runAnalysis(self.y);
		score=0.0;
		for i in range(0,m):
			score=score+abs(maf[i]-self.noisyMAF[i]);
		score=score/self.noise;
		return math.exp(-score);


	def getProbs(self):
		probs=np.asarray([0.0 for i in range(0,self.n)])
		while not self.is_Done():
			score=getScore();
			probs=probs+score*np.asarray(self.y);
			self.next();
		totVal=sum(probs);
		probs=float(self.n1)*probs/totVal;
		return probs;







if __name__=="__main__":
	n=;
	n1=;
	filename=;
	noise=.1
	iter=ItAll(n1,n,filename,snpfile,noise);
	print iter.getProbs();