from Marginals import *;
import numpy as np;
#from matplotlib import pyplot as plt;
from Simple_Tests import gen_Test_Bed as gen;
import math;
from Gen_Fig_Data import *;


##
##THe log of a choose b
##
def logBinom(a,b):
	if b==0 or b==a:
		return 0.0;

	return sum([math.log(i) for i in range(1,a+1)])-sum([math.log(i) for i in range(1,b+1)])-sum([math.log(i) for i in range(1,a-b+1)])

def exactVals(n0,n1,m,lamb):
	vals=[logBinom(n1,k)+logBinom(n0,(n1-k))+m*((k-n1)/float(n1*lamb)-math.log(lamb)) for k in range(0,n1+1)]
	mx=max(vals)
	vals=[i-mx for i in vals]
	prob_n1=sum([k/float(n1)*math.exp(vals[k]) for k in range(0,n1+1)])
	prob_n0=sum([(n1-k)/float(n0)*math.exp(vals[k]) for k in range(0,n1+1)])

	c=n1/float(n1*prob_n1+n0*prob_n0);

	prob_n1=c*prob_n1;
	
	prob_n0=c*prob_n0;
	
	I=[prob_n1 for i in range(0,n1)]
	I.extend([prob_n0 for i in range(0,n0)])

	return I;


def ViewNoise(n0,n1,m,lamb):
	filename="test/test"
	print "Generate"
	gen(filename,n0,n1,m)
	print"Get Exact"
	truth=exactVals(n0,n1,m,lamb)

	print "Estimate Marginal"
	exact=[1.0 for i in range(0,m)]
	marg=Marginals(filename,snpfile="",noise=lamb,err=-1.0,numStep=10000,numSamples=100,params="",analyser="MAF")
	marg.setExact(exact)

	marginals=marg.getMarginals(True);
	
	print marginals
	print truth
	print sum(truth);
	print sum(marginals)

	#print "Plot it!"
	res=[truth[i]-marginals[i] for i in range(0,len(truth))]
	#plt.show()
	return res;



def loadData_hist(filename):
	fil=open(filename)
	line=fil.readline()
	s=line.strip("[")
	s=s.strip("]\n")
	s=s.split(",")
	s=[float(i) for i in s]
	return s

if __name__=="__main__":
	filename="data_for_fig/Accuracy.txt"
	draw=True
	n0=950;
	n1=50;
	m=10;
	lamb=.1
	#print exactVals(100,20,10,.001)
	if not draw:
		res=ViewNoise(n0,n1,m,lamb)
		fil=open(filename,"w");
		fil.write(str(res))
		fil.close()
	else:
		xlab="Error in Estimate";
		ylab="Number of Trials";
		y=loadData_hist(filename)
		drawFig([],y,xlab,ylab,useHist=True)







