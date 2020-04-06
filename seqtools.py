# These scripts offer support for working with sequencing data and annotation information.

from Bio import Seq
from BCBio import GFF

# Tool for loading the entire yeast genome into memory.
def makeGFFlist(GFFgen):
	GFFlist={}
	for chr in GFFgen:
		GFFlist[chr.id]=chr
	return GFFlist

# Table generator for UTR start and end points (like a GFFlist for UTRs).
# utrgffgen is the gffgen for 3'UTR from Nagalakshmi for R64.
def makeutrtable(utrgffgen):
	table={}
	for chr in utrgffgen:
		for feature in chr.features:
			table[feature.id[:-5]]=[feature.location.start.position,feature.location.end.position]	
	return table

# This is a  wrapper for givegene
def givegene(chromosome,feature,GFFs,counts,bp,goodgenes):
	return(givegene_ScEc(chromosome,feature,GFFs,counts,bp,goodgenes))
	
# This is the givegenefunction. It pulls sequence from the GFF and sequencing data from
# the density file.
# It performs splicing 
# It eliminates bad annotations and overlapping and dubious genes (goodgenes=1) 
# goodgenes=2 allows overlapping genes.
# It is called by other functions.
def givegene_ScEc(chromosome,feature,GFFs,counts,bp,goodgenes):
	
	if type(GFFs)==list:		
		GFFlist=GFFs[0]
		utrtable5=GFFs[1]
		utrtable3=GFFs[2]
	else:
		GFFlist=GFFs
		utrtable5={}
		utrtable3={}
		#print "Warning: UTRs of neighboring features not considered in overlap check."
				
	if type(bp)!=int:
		if len(bp)==2:
			bp.append(0)
	else:
		bp=[bp,bp,0]
	
	if bp[0]<0 or bp[1]<0:
		print "Error - bp is negative!"
		return [-1,-1,-1,-1,-1]
	
	if(GFFlist[chromosome].features[feature].strand==1):
		strand=0
	else:
		strand=1
	
	# Eliminate features with problems for analysis.
	if goodgenes > 0:
		if (GFFlist[chromosome].id == 'chrMito' or GFFlist[chromosome].id == '2-micron'):
			return [-1,-1,-1,-1,-1]
		# Excludes features marked as "region" rather than "gene".
		if (GFFlist[chromosome].features[feature].type!='gene'):
			return [-1,-1,-1,-1,-1]
		# Drop dubious genes.
		if (GFFlist[chromosome].features[feature].qualifiers.has_key("orf_classification")):
			if (GFFlist[chromosome].features[feature].qualifiers["orf_classification"][0]=="Dubious"):
				return [-1,-1,-1,-1,-1]
		# These genes excluded because annotation is incorrect.
		badgenelist=['YGL033W','YJR120W']
		for badgene in badgenelist:
			if GFFlist[chromosome].features[feature].id==badgene:
				return [-1,-1,-1,-1,-1]
		
	# Need to get start and end of CDS
	start=GFFlist[chromosome].features[feature].location.start.position
	end=GFFlist[chromosome].features[feature].location.end.position
	
	#check overlap. Overlapping genes going in different directions will not be considered overlapping.
	if goodgenes == 1:
		feats=neighbors(GFFlist,chromosome,feature)
		prevfeat=feats[0]
		nextfeat=feats[1]
		prevend=GFFlist[chromosome].features[prevfeat].location.end.position
		nextstart=GFFlist[chromosome].features[nextfeat].location.start.position
		prevdir=GFFlist[chromosome].features[prevfeat].strand==GFFlist[chromosome].features[feature].strand
		nextdir=GFFlist[chromosome].features[nextfeat].strand==GFFlist[chromosome].features[feature].strand
			
		prevfeatid=GFFlist[chromosome].features[prevfeat].id
		nextfeatid=GFFlist[chromosome].features[nextfeat].id
		if utrtable5.has_key(prevfeatid):
			prev5=utrtable5[prevfeatid][1]-utrtable5[prevfeatid][0]
		else:
			prev5=0
		if utrtable3.has_key(prevfeatid):
			prev3=utrtable3[prevfeatid][1]-utrtable3[prevfeatid][0]
		else:
			prev3=0

		if utrtable5.has_key(nextfeatid):
			next5=utrtable5[nextfeatid][1]-utrtable5[nextfeatid][0]
		else:
			next5=0
		if utrtable3.has_key(nextfeatid):
			next3=utrtable3[nextfeatid][1]-utrtable3[nextfeatid][0]
		else:
			next3=0

		if strand==0:
			if prevdir==True:
				prevend+=prev3
			if nextdir==True:
				nextstart-=next5
		else:
			if prevdir==True:
				prevend+=prev5
			if nextdir==True:
				nextstart-=next3
			
		# We don't know 5' and 3' end here, so need to determine that.
		if strand==0:
			if (start-bp[0])<prevend and prevdir==True:
				return[-2,-2,-2,-2,-2]
			if (end+bp[1])>nextstart and nextdir==True:
				return[-2,-2,-2,-2,-2]
		else:
			if (start-bp[1])<prevend and prevdir==True:
				return[-2,-2,-2,-2,-2]
			if (end+bp[0])>nextstart and nextdir==True:
				return[-2,-2,-2,-2,-2]
	
	splicedseq=Seq.Seq('')
	splicedcounts=[]

	#Assemble gene by splicing out introns here.
	for item in GFFlist[chromosome].features[feature].sub_features:
		if item.type == 'CDS':
			start_feat = int(item.location.start.position)
			end_feat = int(item.location.end.position)
			splicedseq+=(GFFlist[chromosome][start_feat:end_feat]).seq
			splicedcounts+=counts[strand][chromosome][start_feat:end_feat]		
			
	# Check to see if there were no CDS entries:
	if splicedcounts==[] or str(splicedseq)=='':
		return [-1,-1,-1,-1,-1]
	
	# We'll want to swap ends if negative strand:
	if strand==1:
		bp=[bp[1],bp[0],bp[2]]
		
	# Add on extra bit on front and back.
	splicedcounts=counts[strand][chromosome][start-bp[0]:start]+splicedcounts
	splicedseq=GFFlist[chromosome][start-bp[0]:start].seq+splicedseq
	splicedcounts+=counts[strand][chromosome][end:end+bp[1]]
	splicedseq+=GFFlist[chromosome][end:end+bp[1]].seq
	
	# Perform shift for ribosome counts.
	if bp[2]!=0:
		if strand==0:
			if bp[2]>0:
				splicedcounts=counts[strand][chromosome][start-bp[0]-bp[2]:start-bp[0]]+splicedcounts[:-bp[2]]
			else:
				splicedcounts=splicedcounts[-bp[2]:]+counts[strand][chromosome][end+bp[1]:end+bp[1]-bp[2]]
		else:
			if bp[2]>0:		
				splicedcounts=splicedcounts[bp[2]:]+counts[strand][chromosome][end+bp[1]:end+bp[1]+bp[2]]
			else:
				splicedcounts=counts[strand][chromosome][start-bp[0]+bp[2]:start-bp[0]]+splicedcounts[:bp[2]]
	
	# RC for -1 strand.
	if strand==1:
		splicedcounts.reverse()
		splicedseq=splicedseq.reverse_complement()
	
	# Shifting created a problem, probably because the shift was longer than some genes.
	if len(splicedcounts)!=len(splicedseq):
		return [-1,-1,-1,-1,-1]
		
	
	return [splicedcounts,splicedseq,-1,-1,-1]	# Note that items 2-4 in this list are always -1 since they are not used, for future functionalities.

	

# Helper function returns feature numbers of nearest genes that are nondubious.
def neighbors(GFFlist,chromosome,feature):
	if feature<0 or feature >=len(GFFlist[chromosome].features):
		print "Illegal feature given to neigbors program."
		return[-1,-1]
	
	i=1
	checkvar1=True
	checkvar2=True
	while(checkvar1==True or checkvar2==True):
		if((feature-i)<=0 and checkvar1==True):
			checkvar1=False
			prevfeat=0
		if((feature+i)>=len(GFFlist[chromosome].features) and checkvar2==True):
			checkvar2=False
			nextfeat=0
		if(checkvar1==True):
			if GFFlist[chromosome].features[feature-i].qualifiers.has_key("orf_classification"):
				if (GFFlist[chromosome].features[feature-i].qualifiers["orf_classification"][0]!='Dubious'):
					dub=False
				else:
					dub=True
			else:
				dub=False
			if (GFFlist[chromosome].features[feature-i].type=='gene' and dub==False):
				prevfeat=feature-i
				checkvar1=False
		
		if(checkvar2==True):
			if GFFlist[chromosome].features[feature+i].qualifiers.has_key("orf_classification"):
				if (GFFlist[chromosome].features[feature+i].qualifiers["orf_classification"][0]!='Dubious'):
					dub=False
				else:
					dub=True
			else:
				dub=False			
			if (GFFlist[chromosome].features[feature+i].type=='gene' and dub==False):
				nextfeat=feature+i
				checkvar2=False
		i+=1
		
	return [prevfeat,nextfeat]	   

# Give it chrom pos and it tells you the mrna position.
def chrpostomrnapos(chrpos,chrom,featurenum,GFFlist):
	# Make list of features
	sublist=[]
	for subfeature in GFFlist[chrom].features[featurenum].sub_features:
		if subfeature.type=='CDS':
			start=subfeature.location.start.position
			end=subfeature.location.end.position
			sublist.append([start,end])
	prevlength=0
	
	if(GFFlist[chrom].features[featurenum].strand==-1):
		sublist.reverse()
	
	if len(sublist)==0:
		return -1
	
	for item in sublist:
		start=item[0]
		end=item[1]
		length=end-start
		
		if(GFFlist[chrom].features[featurenum].strand==1):
			if chrpos >= start and chrpos < end:		
				mrnapos=prevlength+chrpos-start
				return mrnapos
			else:
				prevlength+=length
		elif(GFFlist[chrom].features[featurenum].strand==-1):
			if chrpos < end and chrpos >= start:
				mrnapos=prevlength+end-chrpos-1		
				return mrnapos
			else:
				prevlength+=length
		else:
			print "Should not be here."
	
	# This is a UTR.
	if(GFFlist[chrom].features[featurenum].strand==1):
		mrnapos=prevlength+chrpos-start-length	# need -length to compensate for prevlength being adjusted.
	else:
		mrnapos=prevlength+end-chrpos-1-length	# need -length to compensate for prevlength being adjusted.
	return mrnapos
