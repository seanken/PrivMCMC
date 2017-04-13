# PrivMCMC

This repository contains the code for the PrivMCMC package. This package takes in genotype and phenotype data, and allows a use to estimate the privacy loss by releasing statistics related to that data. Here, the loss of privacy is measured in terms of posterior probabilities (see our manuscript for details).

Dependencies:

This package has numerous dependencies, namely:
1) pysnptools
2) numpy
3) scipy
4) matplotlib

Running an Analysis with MAF or Log Odds Ratio:

The simplest case for how to use the package is when releasing MAF or Log Odds Ratio. In particular, the package assumes SNP and disease data (present as a binary plink format).

In order to use the package, users should add the directory containing the python files to their path.

A user begins by importing the Marginals package:

import Marginals

The user then creates a marginal object:

marg=Marginals(filename,snpfile,noise,err,numStep=1000,numSamples=1000,params="",analyser="MAF")

In this case, filename is the prefix of the plink file of interest, snpfile is a file with one line for each SNP to be analysed, with the name of that snp on the line. The noise variable is the epsilon noise parameter used in the text, and err is the maximum error. numStep is the number of steps taken in each step of the iteration, numSamples is the number of samples to be drawn. The analyser specifies the type of analysis: either “MAF”, “LogOdds”, “Log” (note this option is still buggy), or a novel analysis (details below). The params parameter is a parameter passed to the analysis—see below.

Having generated the Marginals object, the user can then estimate the marginals:  

results=marg.getMarginals(burnIn=100000)

here, burnIn is the burn-in used by the MCMC algorithm. The returned list, results, has one entry for each individual in the data set. This entry is the estimated posterior probability of that individual being in the case cohort given the released statistics.

Other Statistics:

Our implementation also allows novel statistics, not just MAF and LogOdds. In order to do this, a use must create a Class that extends the Analysis class. This class must implement two methods:

1) setup(self):

The function setups the analysis object (aka does pre computation). 

2) runAnalysis(self,y):

This function takes in a phenotype vector y (which consists of 0’s and 1’s), and outputs the results of the analysis.

Note that the Analysis class has some built in variables, namely self.BED (A bed object from pysnptools.snpreader) containing all genotype information from the bed file, self.y which contains the original phenotype information, self.snps which contains a list of SNPs to be analysed, and self.params, which contains any parameters passed to the constructor (passed using the param parameter).

To see an example of how this is done, look at the MAF_Analysis and Log_Odd_Analysis classes in Analysis.py. More detail will be added at a later time.

Other files:

The remaining files were the ones used to generate the figures used in the manuscript, and some simple test cases/ sanity checks on the different functions in PrivMCMC. We do not give a detailed overview of them here.
