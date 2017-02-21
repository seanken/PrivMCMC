from Analysis import *;
import numpy as np;
import scipy as sp;
import math;
from pysnptools.snpreader import Bed;
from numpy.random import laplace as lap;


class DatasetIterator:
	def __init__(self,filename,snpfile,noise=1.0,err=.01,params="",analyser="MAF",randStart=True,n0=-1,n1=-1):
		if analyser=="MAF":
			self.analyser=MAF_Analysis(filename,snpfile,"",n0,n1);
		if analyser=="Log":
			self.analyser=Log_Analysis(filename,snpfile,params,n0,n1);
		if analyser=="LogOdds":
			self.analyser=Log_Odd_Analysis(filename,snpfile,"",n0,n1)
		self.noise=noise;
		self.err=err;
		self.perm=[i for i in range(0,self.analyser.n)]
		if randStart:
			rand.shuffle(self.perm);
		self.initPerm=[i for i in self.perm]
		temp=self.getTrueValue();
		self.y=self.analyser.getY();

	##
	##Allows user to use non-standard analysis tool
	##
	def setAnalyser(self,analyser):
		self.analyser=analyser;
	
	
	##
	##Returns number people
	##
	def getN(self):
		return len(self.y);
	
	##
	##Restarts
	##
	def reInitialize(self,randStart=True):
		self.perm=[i for i in self.initPerm]
		if randStart:
			rand.shuffle(self.perm);

	def burnIn(self,numSteps=100000):
		for i in range(0,numSteps):
			self.step();


	def noNoise(self):
		exact=self.analyser.getAnalysis();
		print exact

	##
	##Gets value to be released
	##
	def getTrueValue(self):
		exact=self.analyser.getAnalysis();
		m=len(exact);
		for i in range(0,m):
			pert=lap(scale=self.noise);
			while abs(pert)>self.err and self.err>0:
				pert=lap(scale=self.noise);
			exact[i]=exact[i]+pert;
		self.truth=[i for i in exact]
		return [i for i in exact];

	def setExact(self,exact):
		self.truth=exact;
	
	def step(self):
		i=0;
		j=0;
		while self.y[self.perm[i]]==self.y[self.perm[j]]:
			i=rand.randint(0,len(self.y)-1);
			j=rand.randint(0,len(self.y)-1);
		cur=self.analyser.getAnalysis(self.perm);
		permi=self.perm[i];
		self.perm[i]=self.perm[j]
		self.perm[j]=permi;
		try:
			proposed=self.analyser.getAnalysis(self.perm);
		except:
			print("ouch")
			self.perm[j]=self.perm[i];
			self.perm[i]=permi;
			return;
		score1=-sum([abs(cur[k]-self.truth[k]) for k in range(0,len(cur))])/self.noise;
		score2=-sum([abs(proposed[k]-self.truth[k]) for k in range(0,len(cur))])/self.noise;
		if score2>score1:
			return;
		prob=math.exp(score2-score1);
		val=rand.uniform(0,1);
		if val<prob:
			return;
		self.perm[j]=self.perm[i];
		self.perm[i]=permi;






	def sample(self,numSteps=100):
		for i in range(0,numSteps):
			self.step();
		try:
			curVal=self.analyser.getAnalysis(self.perm);
		except:
			self.reInitialize(self,randStart=True)
			return sample(numSteps);
		for i in range(0,len(self.truth)):
			if abs(curVal[i]-self.truth[i])>self.err and self.err>0:
				print "Again!"
				print self.err;
				return self.sample(numSteps);
		return [self.y[i] for i  in self.perm];





def randListSnp(filename,savename):
	print "To be implemented!"










