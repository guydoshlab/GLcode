# Code for characterizing pausing in ribosome profiling data.

import struct
from BCBio import GFF
import listavg
import gentools
import csv
import seqtools
import builddense

# Simple workflow for makeposavg. It averages ribosome profiling data around sites of interest.
# genelist - sites to average, a csv with gene names in col 0, chrom id in 2, feature number in 3, and chromosome position in 4. *** NOTE ***> In 5 is the mRNA position that can be used by commenting that line out below. For 5'UTR it is required.
# GFFfile - path and name of GFF file
# utrgfffile - this is either the GFF for the 5'UTR or 3'UTR, depending on whether your sites are in one of those areas. It can also be 0 if sites are in ORFs.
# seqwin - this is a list of the number of nt upstream and downstream of your feature of interest to include in the average.
# densityfile - This is the path and fileroot to the ribosome profiling data.
# outfilestring - root of the csv file you will have generated
# thresh - generally this is a value for minimal number of rpm counts within the total seqwin to be included. Alternatively, this can signal that the values in the seqwin window should be normalized by the ORF rpm per base. To do this, make thresh a list, with the 2nd value the minimal number of RPKM in the ORF and the first number minimal number rpm in the window.
# riboshift - Amount to shift density before analysis.
def makeposavg_wf0(genelist,GFFfile,utrgfffile,seqwin,densityfile,outfilestring,riboshift,thresh):
	gffgen=GFF.parse(GFFfile)
	GFFlist=seqtools.makeGFFlist(gffgen)
	
	doingUTR5=0
	if type(utrgfffile)==list:
		utrgff=GFF.parse(utrgfffile[0])
		utrtable=seqtools.makeutrtable(utrgff)
		doingUTR5=1
	elif utrgfffile=="0":
		utrtable=0
	else:
		utrgff=GFF.parse(utrgfffile)
		utrtable=seqtools.makeutrtable(utrgff)
	GFFlists=[GFFlist,utrtable,doingUTR5]
	
	counts1p=builddense.readcountsf(densityfile+"_plus_")
	counts1m=builddense.readcountsf(densityfile+"_minus_")
	readcounts=[counts1p,counts1m]

	genesinavg=makeposavg(genelist,GFFlists,seqwin,readcounts,outfilestring,riboshift,thresh)
	print "positions in avg = "+str(genesinavg[0])
	print "positions not in avg because zero count = "+str(genesinavg[1])
	print "positions with zero ORF count = "+str(genesinavg[2])

# This function will take a list with gene names and positions and a window over which to average and then output the average as a binary.
def makeposavg(genelist,GFFlists,seqwin,readcounts,outfilestring,riboshift,thresh):
	zeroORFcount=0
	if type(thresh)==list:
		ORFthresh=thresh[1]
		tempvalthr=thresh[0]
		thresh=tempvalthr
		ORFnorm=1
	else:
		ORFnorm=0
		
	seqcheck=[0,0]
	if seqwin[0]<0:
		seqwin[0]*=-1
		seqcheck[0]=1
	 	
	if seqwin[1]<0:
		seqwin[1]*=-1
		seqcheck[1]=1

	if type(GFFlists)==list:
		GFFlist=GFFlists[0]
		utrtable=GFFlists[1]
		doingUTR5=GFFlists[2]
	else:
		GFFlist=GFFlists
		utrtable=0
		doingUTR5=0
		
	if doingUTR5==1:
		print "Doing 5'UTR. Ensure use of mRNA in code-> mrnaposition=int(pausedict[genename][4]"

	countslist=[]
	f_csv=open(genelist)
	pausedict=listavg.readindict(f_csv)
	avgcounts=[0 for x in range(seqwin[0]+seqwin[1])]
	count=0
	zerocount=0
	for genename in pausedict:
		if genename=="headers":
			continue
		position=pausedict[genename][3]
		chrom=pausedict[genename][1]
		feat_num=pausedict[genename][2]
		
		position=int(position)
		feat_num=int(feat_num)
		
		#######################################################################################
		####### One option should be chosen. See above comments at top of function. ############
		mrnaposition=seqtools.chrpostomrnapos(position,chrom,feat_num,GFFlist)
		#mrnaposition=int(pausedict[genename][4])			

		if utrtable!=0:
			featid=GFFlist[chrom].features[feat_num].id
			if utrtable.has_key(featid):
				UTRlen=utrtable[featid][1]-utrtable[featid][0]
			else:
				UTRlen=0
		else:
			if seqcheck[1]==1:
				UTRlen=seqwin[1]
			else:
				UTRlen=0
			
		if doingUTR5==1:
			UTR5len=UTRlen
			UTRlen=0
		else:
			UTR5len=0
				
		gg=seqtools.givegene(chrom,feat_num,GFFlist,readcounts,[UTR5len,UTRlen,riboshift],1)			
		if gg[0]==-1 or gg[0]==-2:
			continue
		genelength=len(gg[0])
		genecounts=gg[0]
		if utrtable!=0:
			tooclosetostop=(genelength-UTRlen)	# just ORF.
			tooclosetostart=UTR5len
		else:
			tooclosetostop=0
	
		windowscheck=0
		if doingUTR5==1:
			if ((seqcheck[0]==1 or mrnaposition>=seqwin[0]) and (seqcheck[1]==1 or (tooclosetostart-mrnaposition)>=seqwin[1])):
				windowscheck=1
		else:	
			if ((seqcheck[1]==1 or (genelength-mrnaposition)>=seqwin[1]) and (seqcheck[0]==1 or mrnaposition>=(seqwin[0]+tooclosetostop))):
				windowscheck=1
			
		if windowscheck==1:
			loccounts=genecounts[mrnaposition-seqwin[0]:mrnaposition+seqwin[1]]
			ORFcounts=genecounts[UTR5len:genelength-UTRlen]
			
			if sum(loccounts)>thresh:
				count+=1	
			else:
				zerocount+=1
				continue	
				
		else:
			loccounts=[] # I don't think we need this line but doesn't hurt anything.
			continue
		countslist.append([genename,loccounts,ORFcounts])					
	
	if count==0:
		print "No genes to average."
		return [count,zerocount,0]
			
	# Make average
	for gene in countslist:
		genesum=sum(gene[1])
		ORFsum=1000*sum(gene[2])/len(gene[2]) # so in rpkm
		for i in range(len(avgcounts)):
			if ORFnorm==0:
				if genesum!=0:
					avgcounts[i]+=gene[1][i]/float(genesum)		# Equal weighting.
			else:
				if ORFsum>ORFthresh:
					avgcounts[i]+=gene[1][i]/(float(ORFsum)/1000)		# ORF weighting. 
		if ORFnorm!=0:
			if ORFsum<=ORFthresh:
				zeroORFcount+=1
		
		
	for i in range(len(avgcounts)):
		if ORFnorm==0:
			avgcounts[i]/=count
		else:
			avgcounts[i]/=(count+zerocount-zeroORFcount)	
		
	avgs=[avgcounts]

	f2=open(outfilestring,"wb")
	for i in range(0,seqwin[0]+seqwin[1]):
		f2.write(struct.pack("f",float(avgs[0][i])))
	f2.close()
	return [count,zerocount,zeroORFcount]
	

# This script computes a pause score at the input sites of interest by taking the reads locally around the peak and dividing by the average reads in the gene.	
# Simple workflow for makeposstats
# genelist - sites to average, a csv with gene names in col 0, chrom id in 2, feature number in 3, and chromosome position in 4. In 5 is the mRNA position that is also required.
# GFFfile - path and name of GFF file
# utrgfffile - 3'UTR GFF file. Only include this if your sites are in the 3'UTR. Otherwise 0.
# seqwin - this is a list of the number of nt upstream and downstream of your feature of interest to include in reported sequence information.
# pausewin - Number of nt to add to either side (list of 2) of pause peak for computation of pause score. For example, [2,2] would give a peak of 2+2+1 = 5 nt total.
# densityfile - This is the path and fileroot to the ribosome profiling data.
# outfilestring - root of the csv file you will have generated
# riboshift - Amount to shift the ribosome density before computing score. This will vary with length of the motif
# motiflen - This is a variable that is not used in the this version of the script. Set to -1.
def makeposstats_wf(genelist,GFFfile,utrgfffile,seqwin,pausewin,densityfile,outfilestring,riboshift,motiflen):
	gffgen=GFF.parse(GFFfile)
	GFFlist=seqtools.makeGFFlist(gffgen)
	if utrgfffile=="0":
		utrtable=0
	else:
		utrgff=GFF.parse(utrgfffile)
		utrtable=seqtools.makeutrtable(utrgff)
	GFFlists=[GFFlist,utrtable]
	
	counts1p=builddense.readcountsf(densityfile+"_plus_")
	counts1m=builddense.readcountsf(densityfile+"_minus_")
	readcounts=[counts1p,counts1m]
	newdict=makeposstats(genelist,GFFlists,seqwin,pausewin,readcounts,riboshift,motiflen)
	listavg.writedicttoexcel(newdict, outfilestring)
	
def makeposstats(genelist, GFFlists, seqwin, pausewin, readcounts, riboshift,motiflen):
	distfromendofmotif=3	# This is 3 since we are using the start of a codon as point of interest.
	if seqwin[0]%3!=0 or seqwin[1]%3!=0:
		print "Warning, seqwin not multiple of 3. Output translated sequence will be wrong if hits are in frame."
	
	if type(GFFlists)==list:
		GFFlist=GFFlists[0]
		utrtable=GFFlists[1]
	else:
		GFFlist=GFFlists
		utrtable=0

	f_csv=open(genelist)
	pausedict=listavg.readindict(f_csv)
	for genename in pausedict:
		if genename=="headers":
			continue
		position=pausedict[genename][3]
		chrom=pausedict[genename][1]
		feat_num=pausedict[genename][2]
		position=int(position)
		feat_num=int(feat_num)
		mrnaposition=int(pausedict[genename][4])
		
			
		if utrtable!=0:
			featid=GFFlist[chrom].features[feat_num].id
			if utrtable.has_key(featid):
				UTRlen=utrtable[featid][1]-utrtable[featid][0]
			else:
				UTRlen=0
		else:
			UTRlen=0

		gg=seqtools.givegene(chrom,feat_num,GFFlist,readcounts,[0,UTRlen,riboshift],2)			
		if gg[0]==-1:
			continue
		genelength=len(gg[0])
		genecounts=gg[0]
		genesequence=gg[1]
		
		#PRINCIPLE DENOM
		if utrtable!=0:		# Case of looking at a 3'UTR.
			tooclosetostop=(genelength-UTRlen)	# no ORF.
			if UTRlen!=0:
				pause_denom=sum(genecounts[tooclosetostop:genelength])/UTRlen
			else:
				pause_denom=0		
		else:				# Case of looking at an ORF.
			tooclosetostop=0
			pause_denom=sum(genecounts)/genelength
			
		# Get pause score
		if((genelength-mrnaposition)>pausewin[1] and mrnaposition>=(pausewin[0]+tooclosetostop)):
			pause_numerator=sum(genecounts[mrnaposition-pausewin[0]:mrnaposition+pausewin[1]+1]) # Note includes middle base (+1).
			pause_numerator/=(pausewin[1]+pausewin[0]+1)	
			if pause_denom>0:
				pausescore=pause_numerator/pause_denom
			else:
				pausescore=-1
		else:
			pausescore=-1
			pause_numerator=-1
			
		# To convert to rpkm
		pause_numerator*=1000
		pause_denom*=1000	
		
		# Get sequence if available fully inside of window
		if((genelength-mrnaposition)>seqwin[1] and mrnaposition>=(seqwin[0]+tooclosetostop)):
			locseq=genesequence[mrnaposition-seqwin[0]:mrnaposition+seqwin[1]]
			transeq=locseq.translate()
		else:
			locseq="Outsideofgene"
			transeq="Outsideofgene"
			#continue
		
		pausedict[genename].append(locseq)
		pausedict[genename].append(transeq)		# Assumes that seqwin is a factor of 3.
		pausedict[genename].append(pausescore)
		pausedict[genename].append(pause_numerator)	# in rpkm
		pausedict[genename].append(pause_denom)
		
		
		if motiflen==-1:
			continue

	pausedict["headers"].append("LocalSeq")
	pausedict["headers"].append("LocalSeq_trans")
	pausedict["headers"].append("PauseScore")
	pausedict["headers"].append("num_rpkm")
	pausedict["headers"].append("denom_rpkm")

	return pausedict
	



# Simple workflow for motifavg. This function scans the entire transcriptome for all possible sequence motifs of a particular length (nt or amino acid).
# It returns an average plot binary for each and then computes the pause score for each.
# GFFfile - path and name of GFF file
# utr5gfffile - Not used at this time. Input ignored.
# utr3gfffile - Not used at this time. Input ignored.
# counts_filestring - This is the path and fileroot to the ribosome profiling data.
# motifsize - size of motif in nt or amino acids, depending on inframe setting
# inframe - Set to 0 if searching by nt in all frames, set to 1 if searching only main frame in nt, set to 2 if searching main frame by amino acid.
# thresh - This variable is not used at this time. 
# outfilestring - rootfile and path for output of script.
# mismatches - This variable is not used at this time. 
# shift - Amount to shift your ribosome profiling data to align it with the peak of interest.
# windowsize - Number of nt to add to either side of pause peak for computation of pause score. For example, 2 would give a peak of 2+2+1 = 5 nt total.
# avgwindow -  A list of 2, nt distance upstream and downstream of peak of interest to include in the average plot and to use for the denominator of the pause score. Note that motifs positioned at the ORF ends that do not "fit" these windows are excluded.

def motifavg_2_wf(GFFfile,utr5gfffile,utr3gfffile,counts_filestring,motifsize,inframe,thresh,outfilestring,mismatches,shift,windowsize,avgwindow):
	codons={} # Just a placeholder
	counts0=builddense.readcountsf(counts_filestring+"_plus_")		
	counts1=builddense.readcountsf(counts_filestring+"_minus_")
	counts=[counts0,counts1]
	GFFgen=GFF.parse(GFFfile)
	GFFlist=seqtools.makeGFFlist(GFFgen)
	avglist=[]
	outlist=[]
	outlist.append(["motif","na","na","na","na","hitsincluded","na","na","tothits","Pausescore"])
	
	GFFs=GFFlist
	motifdata=motifavg_2_simple(GFFs,motifsize,inframe,thresh,mismatches,codons,shift,counts,windowsize,avgwindow)
		
	for mm in motifdata.keys():
		print mm
		motifdata[mm][1]=0 # Variable not used
		if motifdata[mm][5]>0:
			for i in range(sum(avgwindow)):
				motifdata[mm][0][i]/=float(motifdata[mm][5])	#Normalization taking place.
		motifdata[mm][3]=0	# Variable not used
		
		# Get pause scores:
		numerator=sum(motifdata[mm][0][avgwindow[0]-windowsize:avgwindow[0]+1+windowsize])
		denominator=sum(motifdata[mm][0])
		denominator/=len(motifdata[mm][0])
		numerator/=(2*windowsize+1)
		if denominator!=0:
			pause=numerator/denominator
		else:
			pause=0
		motifdata[mm][10]=pause

		outlist.append([mm]+motifdata[mm][1:9]+[motifdata[mm][10]])
		avglist+=(motifdata[mm][0])	#Concatenate average files.
	
	# WRite out list.
	gentools.writelisttoexcel(outlist,outfilestring)		#Includes new pause scores.
	# Write out avg file.
	favg=open(outfilestring+".bin","wb")
	gentools.writelistbinint(avglist,favg)
	favg.close()

def motifavg_2_simple(GFFs,motif,inframe,thresh,mismatches,codons,shift,countslist,windowsize,avgwindow):
	hardshift=shift	# This is the value that can be used to filter which regions of gene are used. Set to shift for no adjustment. Set to 13 to avoid stop codon regions for 3' end alignment.
	
	motifdata={}	
	shiftdif=shift-hardshift
	if inframe==2:
		motiflen=(int(motif))*3
	elif inframe<=1:
		motiflen=int(motif)
	else:	
		print "illegal input for inframe." 
		return
	GFFlist=GFFs
	print "Motif is: "+str(motif)+"."
	# Call givegene for every gene in genome, one time.	
	for chrom in GFFlist:
		feat_num=0
		print chrom
		for feature in GFFlist[chrom].features:
			
			gg=seqtools.givegene(chrom,feat_num,GFFs,countslist,[0,0,shift],2)	
			if gg[1]==-1 or gg[1]==-2:		# For genes that are considered not allowed because they are dubious, nongene, wrong chrom, etc. 
				feat_num+=1
				continue

			genelen=len(gg[0])
			if genelen==0:
				print "Zero length gene, ERROR!"
			seqlen=genelen
			
			position=0
			while position<(seqlen-motiflen):	
				currentsequence=gg[1][position:position+motiflen]
				if inframe==2:
					currentsequence=str(currentsequence.translate())
				else:
					currentsequence=str(currentsequence)
		
				# Check if motif made for this yet. And increment total count
				if not motifdata.has_key(currentsequence):
					motifdata[currentsequence]=[[0 for x in range(sum(avgwindow))],[],0,[],0,0,0,0,0,[],0,[],[]]
					motifdata[currentsequence][8]+=1
				else:
					motifdata[currentsequence][8]+=1
				
				if(sum(avgwindow)>0):
					# Add in counts for binary output of averaged hits.
						if (-avgwindow[0]+position)>=shiftdif and (-avgwindow[0]+position)>=0 and (avgwindow[1]+position)<=(genelen+shiftdif) and (avgwindow[1]+position)<=genelen:	
							normfactor=sum(gg[0][-avgwindow[0]+position:avgwindow[1]+position])	# So every component of average is weighted equally.
							if normfactor>0:
								for i in range(sum(avgwindow)):
									motifdata[currentsequence][0][i]+=((gg[0][i-avgwindow[0]+position])/normfactor)	#NOTE: output will have the first position of motif (shifted by shift) at the first point after midpoint, ie position 35 if avgwindow ranges 0-69.
								motifdata[currentsequence][5]+=1	# Increment count of averaged hits.
				if inframe>=1:
					position+=3
				if inframe==0:
					position+=1
			feat_num+=1
	return motifdata
