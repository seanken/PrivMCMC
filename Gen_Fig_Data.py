import numpy as np;
import scipy as sp;
from Marginals import *;
from pysnptools.snpreader import Bed;
from Analysis import *;
from matplotlib import pyplot as plt;

##gets either MAF or Logistic analyser for testing
def get_analyser(filename,snpfile,n0,n1,type_analysis,params=""):
	if type_analysis=="MAF":
		return MAF_Analysis(filename,snpfile,"",n0,n1);
	if type_analysis=="Log":
		return Log_Analysis(filename,snpfile,params,n0,n1);

##
##Generates data for all figures except one measuring effect on accuracy due to MCMC
##
def genData(savename,noise_list,bnd_list,snp_lst,numSamp=1000,numStep=1000,burnIn=10000,n0=950,n1=50,type="MAF",param=""):
	fil=open(savename,"w")
	filename="../GWAS/cleanMAF"
	snpfile="data_for_fig/snplist.txt"
	fil.write(str(n0)+" "+str(n1)+" "+type+" "+param+"\n\n");
	for i in range(0,len(noise_list)):
		err=noise_list[i];
		bnd=bnd_list[i];
		numSnp=snp_lst[i];
		fil.write(str(err)+" "+str(bnd)+" "+str(numSnp)+"\n");
		X=Bed(filename);
		snps=[X.sid[i] for i in range(0,numSnp)]
		fil2=open(snpfile,"w");
		for s in snps:
			fil2.write(s+"\n");
		fil2.close()
		analyser=get_analyser(filename,snpfile,n0,n1,type,param)
		marg=Marginals(filename,snpfile,err,bnd,numStep,numSamp,param,type,n0,n1)
		marginals=marg.getMarginals(True);
		for m in marginals:
			fil.write(str(m)+" ");
		fil.write("\n\n");

	fil.close();


##for figures comparing accuracy to privacy
def noiseVSpriv(noise_list,numSnp=25,bnd=-1.0,numSamp=1000,numStep=1000,burnIn=10000,n0=950,n1=50,type="MAF",param=""):
	savename="data_for_fig/Noise_"+str(numSnp)+"_"+str(bnd)+"_"+type+".txt"
	n=len(noise_list);
	bnd_list=[bnd for i in range(0,n)];
	snp_lst=[numSnp for i in range(0,n)]
	genData(savename,noise_list,bnd_list,snp_lst,numSamp,numStep,burnIn,n0,n1,type,param);


##for figures comparing number of statistics released to privacy
def numSnpVpriv(snp_list,noise=.01,bnd=-1.0,numSamp=1000,numStep=1000,burnIn=100000,n0=950,n1=50,type="MAF",param=""):
	savename="data_for_fig/numSnp_"+str(noise)+"_"+str(bnd)+"_"+type+".txt"
	n=len(snp_list);
	bnd_list=[bnd for i in range(0,n)];
	noise_list=[noise for i in range(0,n)]
	genData(savename,noise_list,bnd_list,snp_list,numSamp,numStep,burnIn,n0,n1,type,param);

##Draws figure nicely
def drawFig(x,y,xlab="",ylab="",usebnd=True,useHist=False):
	fig=plt.figure();
	ax=fig.add_subplot(1,1,1);
	if not useHist:
		ax.plot(x,y)
	else:
		ax.hist(y,bins=30)
	ax.spines['top'].set_visible(False);
	ax.spines['right'].set_visible(False);
	#ax.spines['left'].set_visible(False);
	#ax.spines['bottom'].set_visible(False);
	ax.xaxis.set_ticks_position('none');
	ax.yaxis.set_ticks_position('none');
	plt.xlabel(xlab,fontsize=20);
	plt.ylabel(ylab,fontsize=20);
	if usebnd and not useHist:
		plt.ylim([0,1])
	#plt.xlim(xl)
	plt.show()

##loads data saved
def loadData(savename,loc=2):
	fil=open(savename);
	lines=fil.readlines();
	fil.close();
	lines=[l for l in lines if len(l)>5]
	lines=lines[1:]
	split_lines=[l.strip().split() for l in lines]

	long_lines=[s for s in split_lines if len(s)>20];
	

	long_lines=[[float(i) for i in s] for s in long_lines]

	max_lines=[max(s) for s in long_lines]

	num_lines=[float(s[loc]) for s in split_lines if len(s)<20];

	return [num_lines,max_lines];












