import dedup
import os

files=["test.fastq"]

infolder="/Users/username/.../"
outfolder="/Users/username/.../"

# CAUTION - files in outfolder will have the SAME name as the input files.

if not os.path.exists(outfolder):
	os.makedirs(outfolder)
		
for footprintfile in files:
	commandstring = "dedup.writeuniques('"+infolder+footprintfile+"','"+outfolder+footprintfile+"')"

	print ""
	print "Command is:"
	print commandstring
	print ""
	eval(commandstring)	







