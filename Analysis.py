import numpy as np;
import scipy as sp;
import math;
import random as rand;
from statsmodels.discrete.discrete_model import Logit as LR;
import pandas as pd;
import sys;
from statsmodels.tools.tools import add_constant as AC;
from pysnptools.snpreader import Bed;
from pysnptools.snpreader import Pheno;
from numpy.random import laplace as lap;



##
##An abstract class, used by DatasetIterator to find solutions
##
class Analysis:
	##
	##Reads in covariate, pheno and geno info from filename (filename in plink format)
	##
	def __init__(self,filename,snpfile="",params="",n0=-1,n1=-1):
		self.BED=Bed(filename);
		self.pheno=Pheno(filename+".fam");
		self.y=self.pheno.read().val[:,3];
		self.y=self.y-1.0;
		self.params=params;
		n=len(self.y)
		
		if n0>0:
			print "Initiate with n0"
			I0=[i for i in range(0,n) if self.y[i]==0.0]
			I0=I0[:n0]
			I1=[i for i in range(0,n) if self.y[i]==1.0]
			I1=I1[:n1]
			I0.extend(I1);
			self.y=self.y[I0]
			self.BED=self.BED[I0,:]

		try:
			if len(snpfile)>0:
				fil=open(snpfile)
				lines=fil.readlines();
				fil.close();
				self.snps=[l.strip() for l in lines]
			else:
				self.snps=self.BED.sid;
		except:
			print "Error loading SNPs!"
			sys.exit();
		self.setUp();
		self.n=len(self.y)
		print "Number of individuals: "+str(self.n)
		self.Cov=[];
		self.params="";


	##
	##Loads Cov;
	##
	def loadCov(self,covfile):
		print "Not yet implemented!"


	##
	##Set the SNPs!
	##
	def setSNPs(self,snpfile="",SNPs=[]):
		if len(SNPs)>0:
			self.snps=[i for i in SNPs]
		else:
			try:
				fil=open(snpfile)
				lines=fil.readlines();
				fil.close();
				self.snps=[l.strip() for l in lines]
			except:
				print "Error loading SNPs!"
				sys.exit();
		self.setUp();


	##
	##sets up for analysis, specific to each subclass, in this method raises error
	##
	def setUp(self):
		raise NotImplementedError;


	##
	##gets analysis results using permutation
	##
	def getAnalysis(self,perm=[]):
		if len(perm)==0:
			return self.runAnalysis(self.y);
		else:
			return self.runAnalysis([self.y[i] for i in perm])


	##
	##Noisy analysis!
	##
	def getNoisyAnalysis(self,noise,err):
		exact=self.getAnalysis();
		m=len(exact);
		for i in range(0,m):
			pert=lap(scale=noise);
			while pert>err and err>0:
				pert=lap(scale=noise);
			exact[i]=exact[i]+pert;
		return exact;



	##
	##Rns analysis. Implemented in sub classes, assumes analysis is 1-D array
	##
	def runAnalysis(self,y):
		raise NotImplementedError;

	##
	##Returns y
	##
	def getY(self):
		return [i for i in self.y];



##
##For MAF
##
class MAF_Analysis (Analysis):
	##
	##Setups up MAF_Analysis (aka loads genotype data)
	##
	def setUp(self):
		self.X=self.BED[:,self.BED.sid_to_index(self.snps)].read().val;
	
	##
	##Gets MAF of case cohort;
	##
	def runAnalysis(self,y):
		I=[i for i in range(0,len(y)) if y[i]==1]
		maf=np.sum(self.X[I,:],axis=0)/float(2*sum(y))
		return maf;




##
##For Logistic, for the moment without covariates
##
class Log_Analysis (Analysis):
	##
	##Setups up Log_Analysis (aka loads genotype data)
	##
	def setUp(self):
		self.X=self.BED[:,self.BED.sid_to_index(self.snps)].read().val;
		self.m=len(self.snps);
		if len(self.params)==0:
			self.params="Odds";
		#if len(self.Cov)>0:#Add in covariates later!
		self.X=AC(self.X,False);
		print np.sum(self.X,axis=0)
		
	
	
	
	
	
	
	##
	##Gets OR of case cohort;
	##
	def runAnalysis(self,y):
		log_res=[0 for i in range(0,self.m)];
		for i in range(0,self.m):
			I=[i]
			I.extend([-1])
			x=self.X[:,I];
			lr=LR(y,x);
			res_lr=lr.fit(disp=0)
		
			if self.params=="Coef":
				log_res[i]=float(res_lr.params[0]);
			if self.params=="Odds":
				coef=float(res_lr.params[0]);
				log_res[i]=math.exp(coef);
			"""
			if self.params=="pval":
				log_res[i]=;
			if self.params=="logpval":
				pval=;
				if pval>0:
					log_res[i]=-np.log10(pval);
				else:
					log_res[i]=-1.0;
			"""
		
			
		return np.asarray(log_res);







##
#log odd ratio, for the moment without covariates
##
class Log_Odd_Analysis (Analysis):
        ##
        ##Setups up Log_Analysis (aka loads genotype data)
        ##
        def setUp(self):
                self.X=self.BED[:,self.BED.sid_to_index(self.snps)].read().val;
                self.m=len(self.snps);
                #if len(self.Cov)>0:#Add in covariates later!







        ##
        ##Gets OR of case cohort;
        ##
        def runAnalysis(self,y):	
		I=[i for i in range(0,len(y)) if y[i]==1]
		sumAllele_case=np.sum(self.X[I,:],axis=0)
		sumAllele_all=np.sum(self.X,axis=0)
		sumAllele_control=sumAllele_all-sumAllele_case
		n=len(y)
		m=self.m
		ORlist=[]
		
		for i in range(0,m):
			a=sumAllele_case[i]
			b=sumAllele_control[i]
			n1=sum(y);
			n0=n-n1;
			OR=(a/b)/((2*n1-a)/(2*n0-b))
			OR=math.log(OR,2)	
			ORlist.append(OR)
		return ORlist;

