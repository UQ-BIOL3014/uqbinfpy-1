from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import itertools
import operator
"""This python module reads in sam files from RNA-seq experiment and processes them and RNA-seq data"""


def sam_reader(filename):
    """Mandatory fields are QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL
for more info http://samtools.github.io/hts-specs/SAMv1.pdf """ 
    data=[]
    f= open(filename,'r')
    for row in f:
        if row.startswith('@'): # skip the header
            pass
        else:
            info=row.strip().split('\t')
            data.append(info)
    return data 


def base_percentages(reads):
    "reports base percentage  %A,%T,%C,%G "
    all_seqs=[]
    for read in reads:
        seq=read[9]
        seq=[seq[i:i+1] for i in range(0,len(seq),1)]
        for nuc in seq:
            all_seqs.append(nuc)
    counts=dict(Counter(all_seqs))
    nucs=counts.keys()
    freqs={}
    for nuc in nucs:
        freqs[nuc]=float(counts[nuc])/sum(counts.values())
    return freqs
    

def numberofreads(reads):
    """Incremented for every sequence-containing line in the sam file, regardless of whether it represents an alignment.
for some files, this is not actually the number of reads. indeed, this may be a poor name for this stat"""
    return len(reads)


def mapped_reads(reads,paired_end=True):
    """If duplicate tracking was enabled via -D, then this attempts to recapitulate the number of unique, mapped, probe-id's in the original sam file. It is multiplied by 2 for paired-end data with duplicate read id's.
The idea is that if you divide this by the number of reads in the fastq you aligned (possibly from the output of fastq-stats),
you will get an accurate "percentage of reads aligned" statistic.
"mapped" is something with a non-negative position, and a "non-asterisk" cigar string."""
    mapped_reads=[]
    store_reads=[]
    for read in reads:
        if read[3]>0 and read[5]!='*':
            mapped_reads.append(read[0])
            store_reads.append(read)
    mapped=set(mapped_reads)
    list_mapped=list(mapped)
    if paired_end==True:
        mapped=len(mapped)+len(mapped)
    else:
        mapped=len(mapped)
    print "number of mapped reads",mapped
    return store_reads
            
    
def mappedBases(mapped_reads):
    """Total number of mapped bases in sam file"""
    seq=""
    for read in mapped_reads:
        seq=seq+read[9]
    return len(seq)

def forward(mapped_reads):
    """The number of lines in the sam file that were aligned to the "forward" strand. No accounting is done on duplicates."""
    forward=[read for read in mapped_reads if read[9]>0]
    return forward
            

def reverse(mapped_reads):
    """The number of lines in the sam file that were aligned to the "reverse" strand. No accounting is done on duplicates."""
    reverse=[read for read in mapped_reads if read[9]<0]
    return reverse 



########Qualities and STATS


def subgroups(mapped_reads):
    """form groups p<1e-3 one group,1e-3<=p<1e-2 one group,1e-2<=p<1 one group a total of three groups"""
    group1=[]
    group2=[]
    group3=[]
    for read in mapped_reads:
        if int(read[4])>29:
            group1.append(read)
        elif int(read[4])<=29 and int(read[4])>12:
            group2.append(read)
        elif int(read[4])<=12:
            group3.append(read)
        else:
            pass
    print len(group1),"in p<1e-3 group"
    print len(group2),"in 1e-3<=p<1e-2 group"
    print len(group3),"in 1e-2<=p<1 group"
    return group1,group2,group3
               


def dinuc_freq(mapped_reads):
    "reports dinucleotide composition using p(Rho) statistics for overrepresentation"
    all_seqs=[]
    for read in mapped_reads:
        seq=read[9]
        seq=[seq[i:i+1] for i in range(0,len(seq),1)]
        for nuc in seq:
            all_seqs.append(nuc)
    counts=dict(Counter(all_seqs))
    nucs=counts.keys()
    freqs={}
    for nuc in nucs:
        freqs[nuc]=float(counts[nuc])/sum(counts.values())
    all_seqs=[]
    for read in mapped_reads:
        seq=read[9]
        seq=[seq[i:i+2] for i in range(0,len(seq),2)]
        for nuc in seq:
            all_seqs.append(nuc)
    counts=dict(Counter(all_seqs))
    dinucs=counts.keys()
    dinuc_counts={}
    for i in dinucs:
        val=float(counts[i])/sum(counts.values())
        dinuc_counts[i]=val/(freqs[i[0]]*freqs[i[1]]) # p-values
    return dinuc_counts


def PercentReadsAligned(group1,group2,group3,numfastq):
    """Provide a list of mapped_reads and the number of reads in the fastq file"""
    mapped_reads=group1+group2+group3
    Mapped=len(mapped_reads)/float(numfastq)
    Unmapped=1-float(Mapped)
##    print "Mapping stats"
##    print"p<1e-3", len(group1)/float(numfastq)
##    print"1e-3<=p<1e-2",len(group2)/float(numfastq)
##    print "1e-2<=p<1",len(group3)/float(numfastq)
##    print "Unmapped",Unmapped
    labels="p<1e-3","1e-3<=p<1e-2","1e-2<=p<1","Unmapped"
    x=[len(group1)/float(numfastq),len(group2)/float(numfastq),len(group3)/float(numfastq),Unmapped]
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.pie(x,labels=labels,autopct='%1.1f%%', shadow=True)
    plt.title('Mapping stats')
    plt.show()
    return Mapped 
        


def length_stats(group1,group2,group3):
    """returns basic stats relating to the lengths of the reads
Calculations are based on the the length of the (possibly hard-clipped) sequence in the sam file."""
    reads=[group1,group2,group3]
    data=[]
    for i in range(0,len(reads)):
        lengths=[]
        for read in reads[i]:
            if int(read[8])<0:
                length=-1*int(read[8])
            else:
                length=int(read[8])
            lengths.append(length)
        mean_len=np.mean(lengths)
        print "group"+str(i+1)+"mean",mean_len
        max_len=np.max(lengths)
        print "group"+str(i+1)+"max length",max_len
        min_len=np.min(lengths)
        print "group"+str(i+1)+"min length",min_len
        data.append(["group"+str(i+1),mean_len,max_len,min_len])
    return data

def plot_length_distrib(group,name):
    """distribution of lengths of all the sam reads"""
    lengths=[]
    for read in group:
        if int(read[8])<0:
            length=-1*int(read[8])
        else:
            length=int(read[8])
        lengths.append(length)
    ##Visualize length distribution
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    n, bins, patches = plt.hist(lengths,100, normed=0, facecolor='g')
    plt.xlabel("lengths")
    plt.ylabel("number of mapped reads")
    plt.title(name)
    plt.show()

def inv_logit(p):
    return 10**(p/-10) 



def plot_base_composition(reads,sym):
    "reports nucelotide frequencies at each position in the sam sequences"
    #DNA_Alphabet=["A","C","T","G","N"]
    all_nucs=[]
    for read in reads:
        nucs={}#dictionary to store nucleotide data
        seq=read[9]
        for i in range(0,len(seq)):
            nucs[str(i+1)]=seq[i]
        all_nucs.append(nucs)
    all_items=[]
    counts=[]
    pos=range(1,len(seq)+1) 
    for dicts in all_nucs:
        for item in dicts.items():
            all_items.append(item)
    all_items.sort(key=operator.itemgetter(0))
    groups= [map(operator.itemgetter(1),list(group)) for key, group in itertools.groupby(all_items, operator.itemgetter(0))]
    for group in groups:
        counts.append(group.count(sym))
    print counts
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.bar(pos,counts,facecolor='g')
    plt.xlabel("Position")
    plt.ylabel("number of mapped reads")
    plt.title(sym)
    plt.show()       
    return counts 
        
#####################################################
#Transcript reader

def transcript_reader(filename):
    data={}
    f= open(filename,'r')
    for row in f:
        if row.startswith('Transcription'): # skip the header
            pass
        else:
            info=row.strip().split('\t')
            if len(info)>17:
                data[info[7]]=[float(info[8]),float(info[9]),float(info[10]),float(info[11]),float(info[12]),float(info[13]),float(info[14]),float(info[15]),float(info[16]),float(info[17])]
            else:
                data[info[6]]=[float(info[7]),float(info[8]),float(info[9]),float(info[10]),float(info[11]),float(info[12]),float(info[13]),float(info[14]),float(info[15]),float(info[16])]
    return data 
            
###############################
#read edgeR results diff express

def Diff_reader_FDR(filename):
    data={}
    f= open(filename,'r')
    for row in f:
        if row.startswith('logFC'): # skip the header
            pass
        else:
            info=row.strip().split('\t')
            if len(info)>17:
                data[info[0]]=(float(info[1]),float(info[2]),float(info[3])) #keep FDR adjusted-pvalues and LOGFC
            else:
                data[info[0]]=(float(info[1]),float(info[2]),float(info[3]))
    return data 
     
    
###########visualise the adjusted p-values
def plot_FDRpval(data):
    """plot distribution of FDR adjusted p-values"""
    vals=[]
    for i,s,v in data.values():
        vals.append(v)
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    n, bins, patches = plt.hist(vals,50, normed=0, facecolor='g')
    plt.xlabel("Adjusted P-values")
    plt.ylabel("number of genes")
    plt.title("FDR adjust.pval plot")
    plt.show()   

def plotMA(data,cutoff=0):
    """MA Plot of logfold change vs logpcm"""
    logfc=[]
    logpcm=[]
    sig_logfc=[]
    sig_logpcm=[]
    for i,s,v in data.values():
        if v>cutoff:
            logfc.append(i)
            logpcm.append(s)
        else:
            sig_logfc.append(i)
            sig_logpcm.append(s)
    plt.plot(logpcm,logfc,'o')
    plt.plot(sig_logpcm,sig_logfc,'ro')
    plt.xlabel("logpcm")
    plt.ylabel("logfc")
    plt.title("Fold change vs abundance")
    plt.show()
    
#######Test Methods

        
t1=sam_reader("/Users/mikael/workspace/binfpy/BIOL3014/prac_5/t1.sam")

# determine the number of reads 
reads=numberofreads(t1)
print "number of reads",reads

#base composition
base=base_percentages(t1)
print base

#obtain the mapped reads 
mapped_reads=mapped_reads(t1,True)
print mapped_reads[0:5] 

#number of mapped bases
num_bases=mappedBases(mapped_reads)
print "number of mapped bases",num_bases



############################################
#Group the mapped_reads                    #
############################################
###get the range of numbers that comprimise the mapping quality 
nums=[]
for read in mapped_reads:
    nums.append(read[4])
nums=set(nums)
print "get a feel for the range of mapping qualities in this sam file", sorted(nums)
###Get the probability MAPQ=-10*log10*Pr{mapping position is wrong} 
for num in nums:
    score=inv_logit(int(num))
    print "MAPQ and probability"
    print num,score
    
group1,group2,group3=subgroups(mapped_reads)

#dinuc frequency of the mapped reads

nuc=dinuc_freq(group1)
print "dinucleotide frequency of mapped reads(p<1e-3)",nuc 

#get the percentage of reads aligned need to know number of entries in fastq file 
percent=PercentReadsAligned(group1,group2,group3,144126)
print percent


stats=length_stats(group1,group2,group3)
print stats 

#plot the length of all three subgroups
plot_length_distrib(group1,"p<1e-3")
plot_length_distrib(group2,"1e-3<=p<1e-2")
plot_length_distrib(group3,"1e-2<=p<1")

#plot nucleotide composition along the mapped read
data=plot_base_composition(group1,'A')
data=plot_base_composition(group1,'T')
data=plot_base_composition(group1,'C')
data=plot_base_composition(group1,'G')

##read transcripts processed

diff=Diff_reader_FDR("/Users/mikael/workspace/binfpy/BIOL3014/prac_5/Diff_t10.txt")

#plot adjusted p-value by FDR
plot_FDRpval(diff)
#plot MA of diff expression data showing significant points
plotMA(diff,0.005)

print "**************significant gene products at FDR adjust.pval<0.005*********************" 
for k,v in diff.iteritems():
    if v[2]<0.005:
        print k

        



