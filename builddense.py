# This script is used to manage "density" files that other scripts use as input.

import gentools
import struct

# Function to convert wigfiles to a density file.
# This handles 2 types of wigs (fixed and variable). 
# Fixedstep assumed to be stepsize 1 and start value to be 1. 
# chrlen is a dictionary of chromosomes for variablestep wigs. It specifies the length of each chrom. If Fixedstep, ignored.
def wigtocounts(wigfile,chrlen):
	resultlist={}
	f=open(wigfile)
	chrom=""
	steptype=-1			# Fixed is 1, Variable is 0, and not assigned is -1.
	for line in f:		
		if line[0:9]=='fixedStep':	
			steptype=1
			if line[11:17]=="chrom=":
				chrom=gentools.parsenext(line[17:],'  ')
				chrom=cleanchrom(chrom)  # Remove any end of lines or quotations.
				resultlist[chrom]=[]
			else:
				print "Error - no chrom."
				return -1
			
		elif line[0:12]=='variableStep':
			steptype=0
			if line[14:20]=="chrom=":
				chrom=gentools.parsenext(line[20:],'  ')
				chrom=cleanchrom(chrom)	# Remove any end of lines or quotations.
				resultlist[chrom]=[float(0) for x in range(chrlen)]
			else:
				print "Error - no chrom."
				return -1
			
		else:
			if steptype==1:
				resultlist[chrom].append(float(line))
			if steptype==0:
				
				resultlist[chrom][int(gentools.parsenext(line,'	'))-1]=float(gentools.parselast(line,'	'))
			else:
				continue	# This happens on early comment lines.
	f.close()
	return resultlist

# Support function for converting wig to countslist function above.
def cleanchrom(chrom):
	if chrom[-1]=="\n":
		chrom=chrom.strip()
	if chrom[0]=='\"':
		chrom=chrom[1:-1]
	return chrom

# Function to read in counts files. 
def readcountsf(filestring):
    import struct
    keys=[]
    resultlist={}
    f2=open(filestring+"keys","r")
    for line in f2:
        keys.append(line.rstrip('\n'))  # rstrip takes off the end of line character.
    for chrom in keys:
        resultlist[chrom]=[]
        with open(filestring+chrom,"rb") as f:
            nextval = f.read(4) 
            while nextval != "": 
                resultlist[chrom].append(struct.unpack("f",nextval)[0])    # This returns a tuple, not an int, so need to take 1st item.
                nextval = f.read(4)
	f2.close()
    return resultlist
    
# Function to write counts files to disk.
def writecountsf(resultlist,filestring):  #Resultlist  is the list of counts, filestring is the file prefix for each chr to be written.
    import struct
    f2=open(filestring+"keys","w")
    for chrom in resultlist.keys():
        f=open(filestring+chrom,"wb")
        for position in resultlist[chrom]:
            f.write(struct.pack("f",position))
        f.close()   
        f2.write(chrom+"\n")
    f2.close()  