from Gen_Fig_Data import *;
import sys;
import math;


if __name__=="__main__":
	argv=sys.argv;
	if len(argv)==1:
		noiseVSpriv(noise_list=[.01*i for i in range(1,11)],numSnp=25,bnd=-1.0,numSamp=1000,numStep=10000,burnIn=100000,n0=950,n1=50,type="MAF",param="");
	elif argv[1]=="DP":
		savename="data_for_fig/Noise_25_-1.0_MAF.txt"
		xlab="Odds Ratio";
		ylab=r"Level of Differential Privacy ($e^{\epsilon}$)";
		[x,y]=loadData(savename,loc=0)
		y=[math.exp(i*50.0) for i in y]
		x=[(i*1000.0/50.0) for i in x]
		drawFig(x,y,xlab,ylab)
	else:
		savename="data_for_fig/Noise_25_-1.0_MAF.txt"
		xlab="Expected Error";
		ylab="Disclosure Risk";
		[x,y]=loadData(savename,loc=0)
		drawFig(x,y,xlab,ylab)