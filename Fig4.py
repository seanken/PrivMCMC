from Gen_Fig_Data import *;
import sys;


if __name__=="__main__":
	argv=sys.argv;
	if len(argv)==1:
		#numSnpVpriv(snp_list=[1,2,3,10,20,50],noise=1.0,bnd=-1.0,numSamp=100,numStep=1000,burnIn=1000,n0=150,n1=50,type="LogOdds",param="Odds");
		numSnpVpriv(snp_list=[i*10 for i in range(1,6)],noise=.1,bnd=-1.0,numSamp=1000,numStep=10000,burnIn=100000,n0=950,n1=50,type="LogOdds",param="");
	else:
		savename="data_for_fig/numSnp_0.1_-1.0_Log.txt"
		xlab="Number of Values Shared";
		ylab="Disclosure Risk";
		[x,y]=loadData(savename)
		drawFig(x,y,xlab,ylab)
