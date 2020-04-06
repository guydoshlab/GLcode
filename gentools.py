# This is a set of general use functions for python.
import math
import random
import csv

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

