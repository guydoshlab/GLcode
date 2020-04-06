# This file contains code for creating average (metagene) plots and quantitating genes.

from BCBio import GFF
import struct
import builddense
import seqtools
import csv

# This is the workflow wrapper that does the quantitation of the genes for ORF, 5'UTR, and 3'UTR.
# Inputs:
# filebase - path and name of output file base
# counts_filestring - path to input density files and file base
# bp5/3 - Extra distance to add onto annotated distance of UTRs. If ignoreutr is set to 1, these are the UTR lengths.
# ignoreutr5/3 - 1 to ignore UTR annotation, 0 to use it and genes without annotation are not quantitated.
# shift - How much to shift density data.
# filtermodule - A list of 2 elements. The first element is usually a -1. This causes the termination and stacked scores to not be normalized by ORF reads. Put 0 otherwise. The second element should normally be 0. However, if you want to only quantitate certain genes, make 2nd element a path to the list of genes. 
# thresh - Minimal rpkm to be considered in quantitation.
# GFFgen_filename - GFF path and filename
# utrgfffilename/utr5gfffilename - path and filename to GFFs for 5' and 3' UTRs.  

def totalquant_wf(filebase,counts_filestring,bp5,bp3,ignoreutr5,ignoreutr3,shift,filtermodule,thresh,GFFgen_filename,utrgfffilename,utr5gfffilename):

	GFFgen=GFF.parse(GFFgen_filename)
	GFFlist=seqtools.makeGFFlist(GFFgen)

	counts0=builddense.readcountsf(counts_filestring+"_plus_")		
	counts1=builddense.readcountsf(counts_filestring+"_minus_")
	counts=[counts0,counts1]
	
	utrgffgen=GFF.parse(utrgfffilename)
	utrtable=seqtools.makeutrtable(utrgffgen)
	
	utrgffgen=GFF.parse(utr5gfffilename)
	utrtable2=seqtools.makeutrtable(utrgffgen)
	
	mgl=makegenelist(counts,GFFlist,utrtable,utrtable2,bp5,bp3,ignoreutr5,ignoreutr3,shift,filtermodule,thresh)
	writedicttoexcel(mgl,filebase+"_genelist")

# This is the workflow that creates the average (or "metagene") plots.

# Inputs are:
# regionlength5/3 - How many nt to include in average up and downstream of point of interest.
# equalweight - set to 1 if all genes are to be weighted equally in average, otherwise 0.
# alignpos - defines how genes are aligned. 0 means align 5'end of transcript, 1 means start codon, 2 means stop codon, 3 means 3' end of mRNA as defined with bp5, bp3.
# regionlength5 and 3 -  the distances up and down of the alignment position that are considered in the average. Note that UTR inputs need to be -1 or the bp distances (+UTR lengths) big enough to accomodate the regionlengths.
# filebase - path and name of output file base
# counts_filestring - path to input density files and file base
# bp5/3 - Extra distance to add onto annotated distance of UTRs. If ignoreutr is set to 1, these are the UTR lengths.
# ignoreutr5/3 - 1 to ignore UTR annotation, 0 to use it and genes without annotation are not included.
# shift - How much to shift density data in the alignment.
# filtermodule - Put 0 to skip this. To only include certain genes, set this to a path and filename of genes to include. To only include certain stop codons, make this a list of stop codons. To only include certain Kozak contexts, make this a list of Kozak contexts.
# thresh - Minimal rpkm to be considered in the average.
# GFFgen_filename - GFF path and filename
# utrgfffilename/utr5gfffilename - path and filename to GFFs for 5' and 3' UTRs.  
# goodzone - For alignpos 1 or 2, this is the region just after the start or just before the stop that is ignored if doing equal weighting. It eliminates bias due to start/stop peaks.
					
def totalavg_wf(regionlength5,regionlength3,filebase,counts_filestring,bp5,bp3,ignoreutr5,ignoreutr3,shift,filtermodule,thresh,equalweight,GFFgen_filename,utrgfffilename,utr5gfffilename,alignpos,goodzone):
	
	f=open(filebase+"_avg_"+str(alignpos)+"_","wb")

	GFFgen=GFF.parse(GFFgen_filename)
	GFFlist=seqtools.makeGFFlist(GFFgen)
	if filtermodule=='0':
		filtermodule=0
	counts0=builddense.readcountsf(counts_filestring+"_plus_")		
	counts1=builddense.readcountsf(counts_filestring+"_minus_")
	counts=[counts0,counts1]
	
	utrgffgen=GFF.parse(utrgfffilename)
	utrtable=seqtools.makeutrtable(utrgffgen)

	utrgffgen=GFF.parse(utr5gfffilename)
	utrtable2=seqtools.makeutrtable(utrgffgen)
	gene=makeavggene(regionlength5,regionlength3,counts,GFFlist,utrtable,utrtable2,ignoreutr5,ignoreutr3,bp5,bp3,shift,filtermodule,thresh,alignpos,equalweight,goodzone)
	for i in range(0,len(gene)):
		f.write(struct.pack("f",float(gene[i])))
	f.close()
	
# The main quantitation function.
def makegenelist(counts,GFFlist,utrtable,utrtable2,bp5_0,bp3_0,ignoreutr5,ignoreutr3,shift,filtermodule,thresh):
	genelist={}		# The list we will output.
	missedthresh=0
	illegalgenes=0
	genesinlist=0
	ORFdict=0
	
	if type(filtermodule)==list:
		if filtermodule[1]!=0:
			ff=open(filtermodule[1])
			ORFdict=readindict(ff)
		temporaryvar=filtermodule[0]
		filtermodule=temporaryvar
	else:
		print "error in filter module"
		exit()
	
	if bp5_0<0 or bp3_0<0:
		print "error"
		exit()
	
	# Put headers on genelist
	genelist["headers"]=["alias","chromosome","featurenum","UTR3length","rpkm","extra","extraUTR5","Stackscore","Termscore","3'UTR","5'UTR","Note","Sequence","Startscore"]
	
	for chrom in GFFlist:
		feat_num=0
		for feature in GFFlist[chrom].features:
			if ignoreutr3!=1:
				if utrtable.has_key(feature.id):
					bp3=utrtable[feature.id][1]-utrtable[feature.id][0]+bp3_0
					noutr3=False
				else:
					bp3=bp3_0	
					noutr3=True
			else:
				bp3=bp3_0
				noutr3=False
			
			if ignoreutr5!=1:
				if utrtable2.has_key(feature.id):
					bp5=utrtable2[feature.id][1]-utrtable2[feature.id][0]+bp5_0
					noutr5=False
				else:
					bp5=bp5_0	
					noutr5=True
			else:
				bp5=bp5_0
				noutr5=False			
					
			# Import sequence and counts. 
			gg=seqtools.givegene(chrom,feat_num,[GFFlist,utrtable2,utrtable],counts,[bp5,bp3,shift],1)
			genesequence=gg[1]
			genecounts=gg[0]
			
			# Get rid of dubious genes, genes with overlap of others.
			if (genesequence==-1 or genecounts ==-1 or genesequence==-2 or genecounts==-2):
				feat_num+=1
				illegalgenes+=1
				continue
				
			# Define ORF
			start=bp5
			end=len(genecounts)-bp3
			if start==end:
				print "Error, gene length is 0 for gene "+feature.id
				exit()
				
			# Define a special reduced ORF that does not include that start and stop codons and artefacts that are present there.
			start_special=start+15
			end_special=end-15
			
			# To check for any cases with short genes:
			if end_special<0 or end_special<start_special:
				illegalgenes+=1
				feat_num+=1
				continue
			
			extrastart=end
			extraend=len(genecounts)
			
			# Get UTR sequences. If no UTR info available, then it will just use the value of bp_0 for that amount.
			extragenesequence2=genesequence[0:start] # UTR5	
			extragenecounts2=genecounts[0:start]
			extragenecounts=genecounts[extrastart:extraend] #UTR3
			extragenesequence=genesequence[extrastart:extraend]	

			# Compute rpkm for UTRs.
			totalgenereads=float(1000)*sum(genecounts[start_special:end_special])/len(genecounts[start_special:end_special])
			if len(extragenecounts)!=0:
				totalextrareads=float(1000)*sum(extragenecounts)/len(extragenecounts)			
			else:
				totalextrareads=0
			
			if len(extragenecounts2)!=0:
				totalextrareadsUTR5=float(1000)*sum(extragenecounts2)/len(extragenecounts2)
			else:
				totalextrareadsUTR5=0	
				
			# Threshold low genes.
			if totalgenereads<thresh:
				feat_num+=1
				missedthresh+=1
				continue
			
			#Now define extra2, the termination pause score.
			pausehalfregion=2	# This is an arbitrary hardcoded value that changes depending on application. Assumes a shift of 12 or 13 for 80S data.
			totalextrareads2=sum(genecounts[(end-5)-pausehalfregion:(end-5)+pausehalfregion+1]) 

			#Define extra3, the start codon score. 
			pausehalfregion=2	# This is an arbitrary hardcoded value that changes depending on application.
			totalextrareads3=sum(genecounts[(start-0)-pausehalfregion:(start-0)+pausehalfregion+1])
			
			# 6-19-17 Define stack pause score	
			# The -35 number varies depending on application.
			totalextrareads_stackpeak=sum(genecounts[(end-35)-pausehalfregion:(end-35)+pausehalfregion+1])		
			if filtermodule==-1:
				genesumval=1		
			else:
				genesumval=sum(genecounts[start_special:end_special])/len(genecounts[start_special:end_special])	
			
			if genesumval!=0 and totalextrareads2!=0:
				totalextrareads2=totalextrareads2/genesumval
				totalextrareads_stackpeak/=genesumval
			else:
				totalextrareads2=-10
				totalextrareads_stackpeak=-10
			
			# To ensure that only valid UTRs are averaged, we are putting this in.
			if noutr3==True:
				totalextrareads=-10	# No utr in dataset.
				totalextrareads2=-10
			if noutr5==True:
				totalextrareads3=-10	# No utr in dataset.
				totalextrareadsUTR5=-10
			
			# Filter sequences for those of interest 
			if ORFdict!=0 :
				if ORFdict.has_key(feature.id):
					print "Gene included: "+feature.id
				else:
					feat_num+=1
					continue	
				
			if "Alias" in feature.qualifiers:
				alias = feature.qualifiers["Alias"][0]
			else:
				alias = "NA"
						
			if "Note" in feature.qualifiers:
				note = feature.qualifiers["Note"][0]
			else:
				note = "NA"		
		
			# Output chrom, featurenum, counts, extracounts, alias, note.
			genelist[feature.id]=[]
			genelist[feature.id].append(alias)
			genelist[feature.id].append(chrom)
			genelist[feature.id].append(feat_num)	
			if noutr3==False:
				genelist[feature.id].append(bp3-bp3_0)		
			else:
				genelist[feature.id].append(-1)
			genelist[feature.id].append(totalgenereads)	
			genelist[feature.id].append(totalextrareads)
			genelist[feature.id].append(totalextrareadsUTR5)
			genelist[feature.id].append(totalextrareads_stackpeak)
			genelist[feature.id].append(totalextrareads2)	# Term pause score.
			genelist[feature.id].append(str(extragenesequence))
			genelist[feature.id].append(str(extragenesequence2))	# 5' UTR
			genelist[feature.id].append(note)
			genelist[feature.id].append(str(genesequence))
			genelist[feature.id].append(totalextrareads3)	# Start pause score.

			feat_num+=1	
			genesinlist+=1
	print "Genes below threshold = "+str(missedthresh)
	print "Genes dropped by givegene (overlap, undesirable features, etc.) = "+str(illegalgenes)
	print "Genes in list = "+str(genesinlist)
	
	return genelist

# The main metagene function.
def makeavggene(regionlength5,regionlength3,counts,GFFlist,utrtable,utrtable2,ignoreutr5,ignoreutr3,bp5_0,bp3_0,shift,filtermodule,thresh,alignpos,equalweight,goodzone):
	missedthresh=0
	illegalgenes=0
	genesinlist=0
	tooshortlist=0
	averagegene=[0 for x in range(0,(regionlength5+regionlength3))]

	# Set initial values in case not used.
	gz1="notassigned"
	gz2="notassigned"

	if filtermodule!=0 and type(filtermodule)!=list:
		ff=open(filtermodule)
		ORFdict=readindict(ff)

	if bp5_0<0 or bp3_0<0:
		print "Illegal negative values for bp."
		exit()
	
	# First set up whether there is UTR annotation.
	for chrom in GFFlist:
		feat_num=-1
		for feature in GFFlist[chrom].features:
			feat_num+=1
			if ignoreutr3!=1:
				if utrtable.has_key(feature.id):
					bp3=utrtable[feature.id][1]-utrtable[feature.id][0]+bp3_0
					noutr3=False
				else:
					bp3=bp3_0	
					noutr3=True
			else:
				bp3=bp3_0
				noutr3=False
			
			if ignoreutr5!=1:
				if utrtable2.has_key(feature.id):
					bp5=utrtable2[feature.id][1]-utrtable2[feature.id][0]+bp5_0
					noutr5=False
				else:
					bp5=bp5_0	
					noutr5=True
			else:
				bp5=bp5_0	
				noutr5=False	 
			
			# Determine which type of average is being done (0,1,2,or3)
			# Get proper extensions on gene. Check for cases here where UTR will be too short.
			if alignpos==0:		# Align on 5' ends.
				if noutr5==True:
					illegalgenes+=1
					continue
				elif regionlength3>bp5:
					tooshortlist+=1
					continue
				else:	
					bp5=regionlength5+bp5

			elif alignpos==1:	# Align at start codon.
				if noutr5==True:
					illegalgenes+=1
					continue
				elif regionlength5>bp5:
					tooshortlist+=1
					continue
				else:
					bp5=regionlength5
					
			elif alignpos==2:	#align at stop codon.
				if noutr3==True:
					illegalgenes+=1
					continue
				elif regionlength3>bp3:
					tooshortlist+=1
					continue
				else:
					bp3=regionlength3
					
			elif alignpos==3:	# Align at 3' end.
				if noutr3==True:
					illegalgenes+=1
					continue
				elif regionlength5>(bp3):
					tooshortlist+=1
					continue
				else:
					bp3=regionlength3+bp3
				
			# Get counts
			gg=seqtools.givegene(chrom,feat_num,[GFFlist,utrtable2,utrtable],counts,[bp5,bp3,shift],1)
			genecounts=gg[0]
			genesequence=gg[1]
			# Get rid of dubious genes, nongenes, genes with overlap of others.
			if (genesequence==-1 or genecounts ==-1 or genesequence==-2 or genecounts==-2):
				illegalgenes+=1
				continue

			# Define ORF
			start=bp5
			end=len(genecounts)-bp3
			if start==end:
				print "Error, gene length is 0 for gene "+feature.id
				exit()
			
			# Filter sequences for those of interest, i.e. specific sequences, or GO term. Otherwise, will look for dictionary ORFdict.
			if type(filtermodule)==list:
				if filtermodule[0]=="TAA" or filtermodule[0]=="TGA" or filtermodule[0]=="TAG":			
					hasthisstop=0
					for stopcodon in filtermodule:
						if str(genesequence[end-3:end])==stopcodon:
							hasthisstop=1
					if hasthisstop==0:
						continue
				else:
					# This is a start codon w/ context filter:
					hasthisstart=0
					for startcodon in filtermodule:
						if str(genesequence[start-3:start-2])==startcodon:
							hasthisstart=1
							print "Gene included: "+feature.id
					if hasthisstart==0:
						continue		
			elif filtermodule!=0:
				if ORFdict.has_key(feature.id):
					print "Gene included: "+feature.id
				else:
					continue
					
			# Threshold on total gene rpkm: 
			totalgenereads=float(1000)*sum(genecounts[start:end])/len(genecounts[start:end])
			if totalgenereads<thresh:
				missedthresh+=1
				continue
	
			# For each option, add region of interest to the growing average.
			if alignpos==0:
				genesinlist+=1
				countlist=genecounts[0:(regionlength5+regionlength3)]
				
				if equalweight==1:
					totcounts=sum(countlist)
				else:
					totcounts=1
				# Add to growing average.
				if totcounts!=0:	
					for i in range(len(countlist)): 
						averagegene[i]+=countlist[i]/totcounts
				else:
					genesinlist-=1
			
			elif alignpos==1:
				# Check for enough length on 3' end (5' end already checked above):
				if len(genecounts)<(regionlength5+regionlength3):
					tooshortlist+=1
					continue
				else:
					genesinlist+=1
					countlist=genecounts[0:(regionlength5+regionlength3)]
					
					if equalweight==1: # Normalize to a good region that is defined by the region with ORF reads and excluding the good zone.
						if regionlength3>goodzone:
							gz1=goodzone
						else:
							gz1=0 
						totcounts=sum(countlist[(regionlength5+gz1):])
					else:
						gz1=0 # Goodzone not used.
						totcounts=1
					# Add to growing average.
					if totcounts!=0:	
						for i in range(len(countlist)): 
							# Warn for genes with big spikes.
							if countlist[i]/totcounts>100 and i<(regionlength5+gz1):
								print "Potential spike on avg"+str(alignpos)+" gene +"+feature.id
							averagegene[i]+=countlist[i]/totcounts
					else:
						genesinlist-=1
				
			elif alignpos==2:
				# Check for enough length on 5' end (3' end already checked above):
				if len(genecounts)<(regionlength5+regionlength3):
					tooshortlist+=1
					continue
				else:
					genesinlist+=1
					countlist=genecounts[-(regionlength5+regionlength3):]
					
					if equalweight==1: # Normalize to a good region that is defined by the region with ORF reads and excluding the good zone.						
						if regionlength5>goodzone:
							gz2=goodzone
						else:
							gz2=0
						totcounts=sum(countlist[0:(regionlength5-gz2)])
					else:
						gz2=0 # Goodzone not used.
						totcounts=1
					# Add to growing average.
					if totcounts!=0:	
						for i in range(len(countlist)): 
							# Warn for genes with big spikes.
							if countlist[i]/totcounts>100 and (len(countlist)-i)<(regionlength3+gz2):
								print "Potential spike on avg"+str(alignpos)+" gene +"+feature.id
							averagegene[i]+=countlist[i]/totcounts
					else:
						genesinlist-=1
							
			elif alignpos==3:
				genesinlist+=1
				countlist=genecounts[-(regionlength5+regionlength3):]
				
				if equalweight==1:
					totcounts=sum(countlist)
				else:
					totcounts=1
				# Add to growing average.
				if totcounts!=0:	
					for i in range(len(countlist)): 
						averagegene[i]+=countlist[i]/totcounts
				else:
					genesinlist-=1

	for m in range(len(averagegene)):
		if genesinlist!=0:
			averagegene[m]/=genesinlist
		else:
			print "Error, no genes to average at start."		
	
	print "goodzone for avg1="+str(gz1)			
	print "goodzone for avg2="+str(gz2)
	print "Genes under threshold: "+str(missedthresh)
	print "Genes removed for overlap or no UTR or dubious or other: "+str(illegalgenes)
	print "Genes removed because feature(s) too short: "+str(tooshortlist)
	print "Genes in average: "+str(genesinlist)
	return averagegene
	
# A helper function to make a masterdict
def makemasterdict(GFFgen):
	masterdict={}
	masterdict["headers"]=["alias","chromosome","featurenum","Note"]	
	for chr in GFFgen:
		feat_num=0
		for feature in chr.features:
			# In pombe, UTRs don't have feature.ids so we skip (and we'd skip anyway):
			if feature.id=="":
				continue
		
			if "Alias" in feature.qualifiers:
				alias = feature.qualifiers["Alias"][0]
			elif "Name" in feature.qualifiers:	
				alias=feature.qualifiers["Name"][0]
			elif "external_name" in feature.qualifiers:
				alias=feature.qualifiers["external_name"][0]			
			else:
				alias = "NA"
			if "Note" in feature.qualifiers:
				note = feature.qualifiers["Note"][0]
			elif "description" in feature.qualifiers:
				note = feature.qualifiers["description"][0]
			else:
				note = "NA"				
			masterdict[feature.id]=[alias,chr.id,feat_num,note]
			feat_num+=1	
	return masterdict
	
# A helper function to add a dictionary to a master dictionary for all genes.
# appendcollist - a list of the columns in newdict to append in every time.
# namelist - names of those columns
# Dictionary (newdict) must have the keys of feature id.
def combinetomaster(masterdict,newdict,GFFgen,appendcollist,namelist):
	if len(appendcollist)!=len(namelist):
		print "error in masterlist combiner."
	for colname in namelist:
		masterdict["headers"].append(colname)
	for chr in GFFgen:
		for feature in chr.features:
			# In pombe, UTRs don't have feature.ids so we skip (and we'd skip anyway):
			if feature.id=="":
				continue
			for colnum in appendcollist:
				if newdict.has_key(feature.id):
					masterdict[feature.id].append(newdict[feature.id][colnum])
				else:
					masterdict[feature.id].append(-10)
	return masterdict

# Function to write dicts out to csv format.
def writedicttoexcel(genelists,filestring):
	import csv
	writer = csv.writer(open(filestring+".csv", "wb"),delimiter=',')
	if type(genelists)!=list:
		genelists=[genelists]		
	for genelist in genelists:		
		
		if(genelist.has_key("headers")):		
			headerrecord=[]
			headerrecord.append("headers")
			for field in genelist["headers"]:
				headerrecord.append(field)
			writer.writerow(headerrecord)
		for gene in genelist.keys():
			generecord=[]
			generecord.append(gene)
			
			if type(genelist[gene])==list:
				for field in genelist[gene]:
					generecord.append(field)
			else:
				generecord.append(genelist[gene])
			if gene=="headers":				
				continue
			writer.writerow(generecord)

# Make a csv back into a dictionary. 
def readindict(f):
	previousgene=""
	counter=1
	filegen=csv.reader(f,delimiter=',')
	output = {}
	for gene in filegen:
		if gene[0]==previousgene:
			modgenename=gene[0]+"_"+str(counter)	# Note this is assuming they are one after another.
			counter+=1
		else:
			modgenename=gene[0]
			counter=1
		output[modgenename]=[]
		for column in gene[1:]:
			output[modgenename].append(column)
		previousgene=gene[0]

	return output