# This script makes a gene avg for a set of files. Output is a csv file.
import listavg
from datetime import datetime
from BCBio import GFF
import os

# This along with code at end prints out this file to the output so you have it.
thiscode=open("avggene_wf_run.py")
printcode=thiscode.read()
thiscode.close() 

####### DO not change above this line. ##################################################
#########################################################################################
# Input filenames and paths for up to 3 sets of files in different folders.
pathroot3=""
pathroot2=""
pathroot="/Users/username/.../"
files3=[]
files2=[]
files=['example1','example2']

# Input transcriptome sequence and annotations
gfffile="/pathtofiles/yeast.gff"
utr3gfffile="/pathtofiles/yeast_3UTRc.gff"
utr5gfffile="/pathtofiles/yeast_5UTRc.gff"

# Path to output files
outpath="/Users/username/.../"

# Output filenames prefix - what your output file will be called.
outprefix="ex1_"

# To be included, total ORF reads cannot be less than this rpkm threshold. Especially important for equal weighting.
thresh=1

# This shift applies to the thresholding above and goodzone normalization below. 
# The plots themselves will also be shifted but you can shift them back in your plotting program if desired.
# Typically -13 for 5'end alignment. For 3' end, -16.
shift=-15

# This parameter controls which average you do. 0 is for centering at transcript starts, 1 for start codons, 2 for stop codons, and 3 for transcript ends (start of poly(A) tail).
avgregion=2

# These the lengths of the regions to include in the average upstream and downstream of the center position for each of the 4 potential average regions.
# Typically we use 25,50 for region=0, 100,300 for region=1, 300,100 for region=2, and 50,25 for region=3.
regionlength=[300,100]

# Use of UTR annotation:
# Indicate here whether you wish to ignore UTR annotation (and use all genes) or only use only those with UTR annotation. 
# By default, averages at the start/end of the transcript must have annotation so this only applies to start/stop codon averages.
# Typically this is set to 1. Must be 0 for transcript start/end (avgregion=0 or 3).
# For avgregion=0,1 lack of 3'annotation is tolerated. For averegion=2,3 lack of 5' annotation is tolerated.
ignoreutr=1

# This variable controls which genes are eliminated.
# All included = 0. No overlaps, dubious = 1. No dubious, overlaps ok = 2.
# Typically we use 1.
# Note that overlaps are determined by using the annotated UTR information on the neighboring genes, and whatever UTR information is below for the gene of interest.
eliminate=1

# How to weight the average - equal or unequal
# 1 gives equal weight to all genes. 0 weights genes by the level of reads - so a few highly expressed genes heavily dominate.
# For start/stop codon averages, the weighting is done using a "goodzone" of only ORF reads (and tweaked by goodzone parameter below).
equalweight=1

# Region to ignore at end of gene when normalization takes place.
# Only applies if equalweight is 1 and for start/stop codon averages (not transcript start/end averages).
# This does take shift above into account.
# If there is a large spike in ribosome density outside of goodzone, you will get a warning to this effect.
# This can have a big effect on the height of start or stop codon peaks.
# For 80S profiling, 6 is a good value. Larger values, say 15 nt, are appealing but they allow duplicated RP genes to dominate and generate giant start/stop peaks. Use care.
goodzone=6

print "Operation begun at "+str(datetime.now())


####### DO not change below this line. ##################################################
#########################################################################################

goodzone=str(goodzone)
regionlength5=str(regionlength[0])
regionlength3=str(regionlength[1])
ignoreutr=str(ignoreutr)
equalweight=str(equalweight)
avgregion=str(avgregion)
shift=str(shift)
thresh=str(thresh)
eliminate=str(eliminate)
	
filelist=[]
for filename in files:
	filelist.append(pathroot+"/"+filename+"/"+filename)
for filename in files2:
	filelist.append(pathroot2+"/"+filename+"/"+filename)
for filename in files3:
	filelist.append(pathroot3+"/"+filename+"/"+filename)
filelist=str(filelist)

output=outpath+outprefix

print ""
commandstring="listavg.totalavg_wf("+filelist+",'"+output+"','"+gfffile+"','"+utr5gfffile+"','"+utr3gfffile+"',"+regionlength5+","+regionlength3+","+ignoreutr+","+shift+","+equalweight+","+eliminate+","+thresh+","+avgregion+","+goodzone+")"
print "Command is:"
print commandstring
eval(commandstring)	

print "Operation concluded at "+str(datetime.now())
print "Command file is:"
print "*************************************"
print printcode
print
print "*************************************"
print "End command file."
print



