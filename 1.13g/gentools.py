# This is a set of general use functions for python.
import math
import random
import csv
import numpy
import os
import listavg


# This function merges the dictionaries.
# The fillins variable is the number of fillinvals (ie NaNs or -10s) to put in if the 2nd dictionary doesn't have the key. Put -1 to leave these keys out completely.
def mergedictionaries(d1,d2,fillins,fillinval):
	commondict={}
	keys1=d1.keys()
	for key in keys1:
		if d2.has_key(key):
			commondict[key]=d1[key]+d2[key]
		elif fillins!=-1:
			commondict[key]=d1[key]+[fillinval for x in range(fillins)]
	return commondict

# Function to write dict out to csv format.
# Can write out list of many dicts sequentially to the same csv.
def writedicttoexcel(genelists,filestring):
	import csv
	writer = csv.writer(open(filestring+".csv", "wb"),delimiter=',')
	if type(genelists)!=list:
		genelists=[genelists]		### This is also new code -- converting any input dict to a length 1 list.
	for genelist in genelists:		### This is the new code (just this line and changing genelist to genelists in func def.)
		
		if(genelist.has_key("headers")):		# This here presumably to get header as 1st line in csv.
			headerrecord=[]
			headerrecord.append("headers")
			for field in genelist["headers"]:
				headerrecord.append(field)
			writer.writerow(headerrecord)
		for gene in genelist.keys():
			generecord=[]
			generecord.append(gene)
			
			# New Feb 2013, check to see if we have a list or a single value.
			if type(genelist[gene])==list:
				for field in genelist[gene]:
					generecord.append(field)
			else:
				generecord.append(genelist[gene])
			if gene=="headers":				# Skip since we did this above.
				continue
			writer.writerow(generecord)
		
# Make a csv back into a dictionary. f is open file handle. Handles multiple hits for same gene by putting in a _# term for them.
def readindict(f):
	previousgene=""
	counter=1
	filegen=csv.reader(f,delimiter=',')
	output = {}
	for gene in filegen:
		if gene[0]==previousgene:
			modgenename=gene[0]+"_"+str(counter)	# Note this is assuming they are one after another. This is generally true, but may not always be if lists are connected together with cat, for example.
			counter+=1
		else:
			modgenename=gene[0]
			counter=1
		output[modgenename]=[]
		for column in gene[1:]:
			output[modgenename].append(column)
		previousgene=gene[0]

	return output
	
	

#Function to write out any binary list to disk.
def writelistbinint(list,f):    # f is a file handle, i.e. f=open("text.bin","wb")
	import struct
	for element in list:
		f.write(struct.pack("f",element))
	f.close()


# Write out list of lists.
def writelisttoexcel(inlist,filestring):
	import csv
	writer = csv.writer(open(filestring+".csv", "wb"))
	for element in inlist:
		# Check to see if element is not a list.
		if type(element)!=list:
			element=[element]			
		writer.writerow(element)
			
# Note that this doesn't convert strings to integers, etc.
def readexceltolist(filestring):
	f=open(filestring)
	newlist=[]
	reader=csv.reader(f)
	for csvline in reader:
		newlist.append(csvline)
	f.close()
	return newlist
	


		
# Returns all text from the start of the input string up to the occurence of the delimiter.
def parsenext(string,delimiter):
	delimlen=len(delimiter)
	i=0
	strlen=len(string)
	while(i<strlen):
		if string[i:i+delimlen]==delimiter:
			break
		i+=1
		
	return string[0:i]
	
# Same function, now returns last value.
def parselast(string,delimiter):
	delimlen=len(delimiter)
	i=0
	strlen=len(string)
	while(i<strlen):
		if string[i:i+delimlen]==delimiter:
			break
		i+=1
		
	return string[i+delimlen:]


# This function returns a list of all positions for searched element in a list.
def getindexpositions(inlist,element):
	inlist=str(inlist)
	lastpos=-1
	positions=[]
	inlistlen=len(inlist)
	while(lastpos<inlistlen):
		if inlist[lastpos+1:].count(element)!=0:
			curpos=inlist[lastpos+1:].index(element)
		else:
			return positions
		positions.append(curpos+lastpos+1)
		lastpos=curpos+lastpos+1
	return positions
	


# Take list of column lists and use to combine those columns from list of csv docs into a master. Also title lists as well for each column.
#filelist is a list. columnlist and titlelist are lists of lists (even if 1 element).
def combinecsvs(filelist,columnlist,titlelist,outfilestring):
	writerfile=open(outfilestring+".csv", "wb")
	writer = csv.writer(writerfile,delimiter=',')
	
	if len(filelist)!=len(columnlist) or len(filelist)!=len(titlelist):
		print "ERR"
		exit()
	numfiles=len(filelist)
	numcolumns=0
	# Get num of columns:
	for col in columnlist:
		numcolumns+=len(col)
	newcsvcols=[]
	
	for i in range(numfiles):
		if len(columnlist[i])!=len(titlelist[i]):
				print "ERR2"
				exit()
		for j in range(len(columnlist[i])):
			f=open(filelist[i])
			print filelist[i]
			rrr=csv.reader(f)
			nc=extractcolumn(rrr,columnlist[i][j])
			tl=[titlelist[i][j]]
			newcsvcols.append(tl+nc)
			f.close()
	
	# Convert columns to rows.
	for i in range(len(newcsvcols[0])):		# OBviously, we need all columns to be same length.
		newrow=[]
		for j in range(numcolumns):
			newrow.append(newcsvcols[j][i])
		#Write row out.
		writer.writerow(newrow)
	writerfile.close()


# Transpose a csv file and delete original.
def transposecsv(csvfile):
	writerfile=open(csvfile+"_transposed.csv", "wb")
	writer = csv.writer(writerfile,delimiter=',')
	f=open(csvfile+".csv")
	readercsv=csv.reader(f)
	maxrow=0
	for row in readercsv:
		if len(row)>maxrow:
			maxrow=len(row)
	
	columnnum=maxrow
	for column in range(columnnum):
		f.close() 
		f=open(csvfile+".csv")
		readercsv=csv.reader(f)
		col=extractcolumn(readercsv,column)
		writer.writerow(col)
	f.close()
	writerfile.close()
	os.remove(csvfile+".csv")
	os.rename(csvfile+"_transposed.csv",csvfile+".csv")
 
# This function will return the number of mismatches between 2 strings.
# No checking of length for speed.
def mismatchcount(str1,str2):
	mismatches=0
	for (a,b) in zip(str1,str2):
		if a!=b:
			mismatches+=1
	return mismatches
	
