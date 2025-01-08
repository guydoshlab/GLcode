print ""
print "Command file is:"
print "*************************************"
thiscode=open("makeposstats_multidensity_wf_run.py")
printcode=thiscode.read()
thiscode.close()
print
print "*************************************"
print "End command file."
print

import os
import searchcode
from datetime import datetime
import gentools

# This is a simple script for computing the "pause score" for one or more datasets at a list of genome positions that you provide.
# The script is rather experimental so it's best to check the underlying code to familiarize yourself.
# Basically, this computes a "numerator" and "denominator" at the positions of interest you provide in the input. Then it take the ratio. Reported as a csv.
# To compute numerator, it takes the value at the position of interest plus reads on either side of it defined by the pausewin variable. 
# The denominator can be computed many different ways and several options are in the searchcode code comments. It's probably best to try several of them. Default is average counts in the region of interest (i.e. main ORF or UTR).
# The seqwin value determines how much sequence about the site of interest is reported.
# motiflen is a variable that is used for some special features that compute reading frame around the site of interest. See underlying code. We generally set this to -1 if not doing this.

######### Do not change above line.
gfffile="/pathtofiles/yeast.gff"
utr3gfffile="/pathtofiles/yeast_3UTRc.gff"
utr5gfffile="/pathtofiles/yeast_5UTRc.gff"

densityfiles=['exfile1','exfile2']

riboshift="-18"

# put 0 for siteshift if you want to keep it at the 5' end of the motif. Put in motif length if you want 3' end. And subtract 3 if you want 3 in (ie A site since shift shifts for P site).
siteshift="0"

# pausewin is a list of 2 numbers for nt upstream/downstream of the peak to include in the numerator of the pause score.
pausewin=str([0,0])			

# These 2 values define the region to include in the sequence reported in the output csv.
seqwin=str([50,50])

# 5'UTR = 0, ORF = 1, 3'UTR = 2. -2 is special for preAug2020 3'UTR lists that used ORF coordinates.
region="1"

# Where your output csv will go and prefix for the filename.
outfileroot="/Users/username/.../outputs_"

# Sites where you are computing pauses.
genelist="polyK3.csv"


####### Do not change below line.

densityfiles=str(densityfiles)
commandstring="searchcode.makeposstats_wf("+densityfiles+",'"+genelist+"','"+GFFfile+"','"+utrgfffile5+"','"+utrgfffile3+"',"+seqwin+","+pausewin+",'"+outfileroot+"',"+riboshift+","+siteshift+","+region+")"
				
print "Command is:"
print commandstring
eval(commandstring)	
print "code used to run:"
print printcode
print "Operation concluded at "+str(datetime.now())

