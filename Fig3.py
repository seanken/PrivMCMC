from Gen_Fig_Data import *;
import sys;
import math

if __name__=="__main__":
	argv=sys.argv;
	if len(argv)==1:
		numSnpVpriv(snp_list=[i*10 for i in range(1,6)],noise=.01,bnd=-1.0,numSamp=1000,numStep=1000,n0=950,n1=50,type="MAF",param="");
	else:
		snp_list=[i*10 for i in range(1,6)]
		noise=.01
		n1=50;
		n0=950
		n=n0+n1;
		savename="data_for_fig/numSnp_0.01_-1.0_MAF.txt"
		xlab=r"Differential Privacy Parameter ($\epsilon$)";
		ylab="Log Odd Ratio";
		[x,y]=loadData(savename)
		eps=[i*10*.01*n1 for i in range(1,6)]
		y=[math.log(i*n/float(n1)) for i in y]
		drawFig(eps,y,xlab,ylab,usebnd=False)