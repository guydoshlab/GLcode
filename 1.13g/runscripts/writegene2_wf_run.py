import seqtools
import os
from datetime import datetime

thiscode=open("writegene2_wf_run.py")
thiscodetext=thiscode.read()
thiscode.close()

# This code an option for generating write genes of many different input files with varied extension for footprint length.

############## DO NOT CHANGE ABOVE THIS LINE.################

# Use this code normally to list your density files.
densityfiles=str(["/Users/username/.../example1","/Users/username/.../example2"])

# What genes you want:
feature=str(["gene1","gene2"])		# Formal or informal names are okay.

gfffile="/pathtofiles/yeast.gff"
utr3gfffile="/pathtofiles/yeast_3UTRc.gff"
utr5gfffile="/pathtofiles/yeast_5UTRc.gff"	

outpath="/Users/username/.../"

# Here is how much UTR to put on each gene. We used to call these variables "shift" values.
bp5=100
bp3=100

# If this value is set to 1, then the UTRs will be the exact lengths given in the bp values. 
# If set to 0, then the bp values will extend the annotated UTRs. 
# If you don't want UTRs at all, put -1.
manualutrs=1

# Options for intron:
#intron=0; standard way with spliced transcript returned. (Sense direction, 5' to 3')
#intron=1; as in old version of this script, unspliced transcript is returned. (Sense direction, 5' to 3')
#intron=2; unspliced transcript returned in 5' to 3' direction. Unspliced transcript returned on antisense strand in 5' to 3' direction.
#intron=3; same as intron = 2 except the orientation for both sense and antisense will be as they appear on the chromosome. Nothing is flipped. So one will be 5' to 3' and the other 3' to 5'.
intron=0		

############## DO NOT CHANGE BELOW THIS LINE.##################
intron = str(intron)

if manualutrs==1:
	if bp5==0 or bp3==0:
		print "Use manualutrs=-1 to to get no UTR extension. Exiting."
		exit()

if manualutrs==1:
	bp5=str(-bp5)
	bp3=str(-bp3)
elif manualutrs==0:
	bp5=str(bp5)
	bp3=str(bp3)
elif manualutrs==-1:
	bp5="'none'"
	bp3="'none'"
else:
	print "Error in manualutrs"
	exit()
	
print "Operation begun at "+str(datetime.now())

commandstring="seqtools.writegene2_wf("+bp5+","+bp3+","+intron+","+densityfiles+","+feature+",'"+gfffile+"','"+utrgfffile5+"','"+utrgfffile3+"','"+outfilepath+"')"
print "Command is:"
print commandstring
eval(commandstring)	

file=open(outfilepath+"_output.txt","w")
file.write("Command is:\r\n")
file.write(commandstring)
file.write("\r\n\r\n")
file.write("Code in script:\r\n")
file.write(thiscodetext)


print "Operation concluded at "+str(datetime.now())
