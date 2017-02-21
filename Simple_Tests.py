from Marginals import *;
from DataIterator import *;
from Analysis import *;
import numpy as np;
from matplotlib import pyplot as plt
from pysnptools.snpreader import SnpData
from pysnptools.snpreader import Bed


##Generates BED file used in testing
def gen_Test_Bed(filename,n0,n1,m):
	n=n0+n1;
	iid=[["fam_"+str(i),"iid_"+str(i)] for i in range(0,n)]
	sid=["snp_"+str(i) for i in range(0,m)]
	X=[[2.0 for i in range(0,m)] for i in range(0,n1)]
	X.extend([[0.0 for i in range(0,m)] for i in range(0,n0)])
	dat=SnpData(iid=iid,sid=sid,val=X);
	Bed.write(filename,dat);
	fil=open(filename+".fam")
	lines=fil.readlines();
	fil.close();
	fil=open(filename+".fam","w")
	for i in range(0,len(lines)):
		l=lines[i]
		s=l.strip().split();
		if i<n1:
			s[5]="2";
		else:
			s[5]="1";
		l=" ".join(s)+"\n";
		fil.write(l);
	fil.close();





##gets either MAF or Logistic analyser for testing
def get_analyser(filename,snpfile,n0,n1,type_analysis="MAF",params=""):
	if type_analysis=="MAF":
		return MAF_Analysis(filename,snpfile,"",n0,n1);
	if type_analysis=="Log":
		return Log_Analysis(filename,snpfile,params,n0,n1);


##
##Some basic unit tests for Analysis
##
def testAnalysis(filename="test/test",n0=10,n1=100,m=10):
	for type_analysis in ["MAF"]:
		gen_Test_Bed(filename,n0,n1,m);
		##firsts, checks that initialization works
		snpfile="";

		n=n0+n1;
		analyser=get_analyser(filename,snpfile,n0,n1);
		if not analyser.X.shape==(n0+n1,m):
			print "Did not get the correct number of individuals/SNPs!"
			exit(0)
		else:
			print "Correct number!"

		##Checks correct
		exact=analyser.getAnalysis();
		perm=[n-i-1 for i in range(0,m)]
		perm_exact=analyser.getAnalysis(perm);
		if [i for i in exact]!=[1.0 for i in exact] or [i for i in perm_exact]!=[0.0 for i in perm_exact]:
			print "Error in analysis!"
			exit(0)
		else:
			print "Analysis is correct!"

		if not sum(analyser.y)==n1:
			print "Disease Status not read in correctly!"
			exit()
		else:
			print "Disease status passed sanity check"


		print "Passed the basic tests!"



##
##Gets the estimates of marginals
##
def getEst(n0,n1,err):
	print "Not yet implemented!"



##
##Some basic tests for marginals
##
def testMarginals(filename="test/test",n0=10,n1=100,m=10):
	err=.1;
	gen_Test_Bed(filename,n0,n1,m);


	##NEED TO FINISH!


	est=getEst(n0,n1,m,err);
	plt.hist([est[i]-marginals[i] for i in range(0,len(est))],bins=30);
	plt.show();



##
##Some basic tests for DataIterator
##
def testIterator(filename="test/test",n0=10,n1=100,m=10):
	gen_Test_Bed(filename,n0,n1,m);
	try:
		dat_iter1=DatasetIterator(filename,snpfile="",noise=.001,err=.01,params="",analyser="MAF",randStart=True,n0=n0,n1=n1)
		dat_iter2=DatasetIterator(filename,snpfile="",noise=.05,err=.1,params="",analyser="MAF",randStart=True,n0=n0,n1=n1)
		dat_iter3=DatasetIterator(filename,snpfile="",noise=.01,err=-1.0,params="",analyser="MAF",randStart=True,n0=n0,n1=n1)
	except:
		print "Error Initializing DataIterator!"
		exit()

	print "Initialization Worked!!"
	dat_iter1.getTrueValue();
	exact=dat_iter1.truth;
	if sum([abs(e-1.0) for e in exact])/float(len(exact))>.01:
		print "Error in estimating value!"
		exit()
	print exact;
	print "Estimate of truth seems reasonable"


	for i in range(0,100):
		dat_iter2.getTrueValue();

		exact=dat_iter2.truth;
		if abs(exact[0]-1.0)>.1:
			print exact[0];
			print dat_iter2.noNoise()
			print dat_iter2.err
			print "Error with bounded Lap noise!"
			exit()
	print "Bounded Laplace noise seems bounded"

	tot=0.0
	for i in range(0,100):
		dat_iter3.getTrueValue();
		
		exact=dat_iter3.truth;

		tot=tot+sum(np.abs(np.asarray(exact)-1.0))/float(m);
	tot=tot/100.0;

	if abs(tot-.01)>.01:
		print "Noise not as should be!"
		exit()

	print "Noise distribution seems ok"


	try:
		samp=dat_iter2.sample(20)
	except:
		print "Sampling is broken!"


	print "Sampling is ok!!"

	print "Passed all tests!"





if __name__=="__main__":
	testAnalysis()
	testIterator()









