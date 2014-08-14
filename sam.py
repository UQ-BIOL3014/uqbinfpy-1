from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import itertools
import operator
import math
import re
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
        elif int(read[4])<=29 and int(read[4])>17:
            group2.append(read)
        elif int(read[4])<=17:
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

def raw_count_reader(filename):
    data={}
    f= open(filename,'r')
    for row in f:
        if row.startswith('t1'): # skip the header
            pass
        else:
            info=row.strip().split('\t')
            data[info[0]]=[int(info[1]),int(info[2]),int(info[3]),int(info[4]),float(info[5])] #t1,rept1,t10,rept10,len 
    return data 


#####Normalisation methods 


def get_RPKM(data,num_map1,num_map2,num_map3,num_map4):
    """provide number of mapped reads for the two groups of interest and raw count data .This method provides length normalisation to prevent length and total count bias"""
    counts1=[];counts2=[];counts3=[];counts4=[];lengths=[]
    for i,s,ii,ss,v in data.values():
        counts1.append(i)
        counts2.append(s)
        counts3.append(ii)
        counts4.append(ss)
        lengths.append(v)
    rpkms=[];rpkms2=[];rpkms3=[];rpkms4=[];final={}
    #perform RPKM calc
    for i in range(0,len(counts1)):
        if counts1[i]==0:
            rpkm=0
            rpkms.append(rpkm)
        else:
            rpkm=float(counts1[i])/(lengths[i]*(float(num_map1)/10**6))
            rpkms.append(rpkm)
    for i in range(0,len(counts2)):
        if counts2[i]==0:
            rpkm=0
            rpkms2.append(rpkm)
        else:
            rpkm=float(counts2[i])/(lengths[i]*(float(num_map2)/10**6))
            rpkms2.append(rpkm)
    for i in range(0,len(counts3)):
        if counts3[i]==0:
            rpkm=0
            rpkms3.append(rpkm)
        else:
            rpkm=float(counts3[i])/(lengths[i]*(float(num_map3)/10**6))
            rpkms3.append(rpkm)
    for i in range(0,len(counts4)):
        if counts4[i]==0:
            rpkm=0
            rpkms4.append(rpkm)
        else:
            rpkm=float(counts4[i])/(lengths[i]*(float(num_map4)/10**6))
            rpkms4.append(rpkm)
    #return gene names and rpkms 
    for i in range(0,len(data.keys())):
        final[data.keys()[i]]=[float(rpkms[i]),float(rpkms2[i]),float(rpkms3[i]),float(rpkms4[i])]
    return final 

def write_RPKM_data(RPKM_data,filename):
    f=open(filename,'w')
    for i in range(0,len(RPKM_data)):
        f.write("%s\t%d\t%d\t%d\t%d\n"%(RPKM_data.keys()[i],int(RPKM_data.values()[i][0]),int(RPKM_data.values()[i][1]),int(RPKM_data.values()[i][2]),int(RPKM_data.values()[i][3])))
    f.close() 
        



###############Visualize replicates to determine degree of biological variation

def pearson_def(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = np.mean(x)
    avg_y = np.mean(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    return diffprod / math.sqrt(xdiff2 * ydiff2)


def plotreprpkm(rpkm_data,timepoint):
    """plot showing level of agreement between technical replicates for RPKM between replicates and plots coefficient of determination"""
    one=[]
    two=[]
    if timepoint=="t1":
        for i in range(0,len(rpkm_data.values())):
            one.append(int(rpkm_data.values()[i][0]))
            two.append(int(rpkm_data.values()[i][1]))
    else:
        for i in range(0,len(rpkm_data.values())):
            one.append(int(rpkm_data.values()[i][2]))
            two.append(int(rpkm_data.values()[i][3]))
    plt.plot(one,two,'o')
    pcc=pearson_def(one,two)
    R2=pcc**2
    name="""Technical Replicates
R2="""+str(R2)
    m,b= np.polyfit(one,two,1)
    plt.figure(1, figsize=(8,8))
    plt.plot(one, np.array(one)*m +b,'r-') 
    plt.text(3000, max(two)-1000,name , fontsize=12)
    plt.xlabel("RPKM replicate 1")
    plt.ylabel("RPKM replicate 2")
    plt.title(timepoint)
    plt.show()


def plotMAreprpkm(rpkm_data,timepoint):
    """MA Plot of log(RPKM) vs Average log(RPKM) of replicates"""
    m=[]
    a=[]
    if timepoint=="t1":
        for i in range(0,len(rpkm_data.values())):
            y=np.log2(rpkm_data.values()[i][0]+1)-np.log2(rpkm_data.values()[i][1]+1)
            x=(np.log2(rpkm_data.values()[i][0]+1)+np.log2(rpkm_data.values()[i][1]+1))/2
            m.append(y)
            a.append(x)
    else:
        for i in range(0,len(rpkm_data.values())):
            y=np.log2(rpkm_data.values()[i][2]+1)-np.log2(rpkm_data.values()[i][3]+1)
            x=(np.log2(rpkm_data.values()[i][2]+1)+np.log2(rpkm_data.values()[i][3]+1))/2
            m.append(y)
            a.append(x)    
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.plot(a,m,'o')
    plt.axhline(np.mean(m)+1.96*np.std(m),color="green",label="avg diff +1.96(std diff)")
    plt.axhline(np.mean(m)-1.96*np.std(m),color="green",label="avg diff -1.96(std diff)")
    plt.xlabel("Average log(RPKM) of replicates")
    plt.ylabel("Difference in log(RPKM) of replicates")
    plt.legend(loc="lower right")
    plt.title(timepoint)
    plt.show()



def get_cv(data1,condition):
    cvs=[]
    if condition=="t1":
        for i in range(0,len(data1.values())):
            mean = np.mean([data1.values()[i][0],data1.values()[i][1]])
            std=np.std([data1.values()[i][0],data1.values()[i][1]])
            if mean==0.0 and std==0.0:
                pass
            else:
                cv=float(mean+1)/(std+1)
                cvs.append(cv) 
    else:
        for i in range(0,len(data1.values())):
            mean = np.mean([data1.values()[i][3],data1.values()[i][4]])
            std=np.std([data1.values()[i][3],data1.values()[i][4]])
            if mean==0.0 and std==0.0:
                pass
            else:
                cv=mean+1/std+1
                cvs.append(cv)          
    return cvs 
        
        
def get_boxplots(norm,original):
    """distribution of the coeficient of variation across samples (replicates) normalised using the methods provided"""
    bp=plt.boxplot([norm,original],notch=False, patch_artist=True)
    for box in bp['boxes']:
        box.set(color="red")
        box.set(color="blue")
    plt.ylabel("coefficient of variation")
    plt.xlabel("Methods")
    my_xticks = ['RPKM','raw counts']
    x=[1,2]
    plt.xticks(x,my_xticks)
    plt.ylim(0,400)
    plt.show()       
                   
               
def plotavg_cv(norm,original):
    """distribution of the coeficient of variation across samples (replicates) normalised using the methods provided"""
    x=[1,2]
    y=[np.mean(norm),np.mean(original)]
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.bar(x[0],y[0],color="red",label="RPKM")
    plt.bar(x[1],y[1],color="blue",label="Raw counts")
    plt.ylabel("Average coefficient of variation")
    plt.xlabel("Methods")
    ax.xaxis.set_ticklabels([])
    plt.legend(loc="upper right")
    plt.show()


def plotMA(rpkm_data,cutoff=[-1.5,1.5]):
    logfc=[]
    avg_rpkm=[]
    sig_logfc=[]
    sig_avg_rpkm=[]
    logfc2=[]
    avg_rpkm2=[]
    sig_logfc2=[]
    sig_avg_rpkm2=[]
    for i,ii,s,ss in rpkm_data.values():
        fc=np.log2(float(s+1)/(i+1))
        if fc<cutoff[0] or fc>cutoff[1]:
            sig_logfc.append(fc)
            sig_avg_rpkm.append(np.log2(s+1)+np.log2(i+1)/2)
        else:
            logfc.append(fc)
            avg_rpkm.append(np.log2(s+1)+np.log2(i+1)/2)
    for i,ii,s,ss in rpkm_data.values():
        fc2=np.log2(float(ss+1)/(ii+1))
        if fc2<cutoff[0] or fc2>cutoff[1]:
            sig_logfc2.append(fc2)
            sig_avg_rpkm2.append(np.log2(ss+1)+np.log2(ii+1)/2)
        else:
            logfc2.append(fc2)
            avg_rpkm2.append(np.log2(ss+1)+np.log2(ii+1)/2)
    plt.figure(1, figsize=(8,8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.plot(avg_rpkm,logfc,'o',color="blue",label="rep1")
    plt.plot(avg_rpkm2,logfc2,'x',color="blue",label="rep2")
    plt.plot(sig_avg_rpkm,sig_logfc,'o',color="red",label="sig rep1")
    plt.plot(sig_avg_rpkm2,sig_logfc2,'x',color="red",label="sig rep2")
    plt.axhline(cutoff[0],color="orange")
    plt.axhline(cutoff[1],color="orange")
    plt.ylabel("Fold Change (log2)")
    plt.xlabel("Average RPKM (log2)")
    plt.title("MA plot")
    plt.legend(loc="upper left")
    plt.show()
        
   
#######Test Methods

        
t1=sam_reader("/Users/samirlal/Desktop/sam/t1.sam")

# determine the number of reads 
reads=numberofreads(t1)
print "number of reads",reads

#base composition
base=base_percentages(t1)
print base

#obtain the mapped reads 
mapped_read=mapped_reads(t1,True)
print mapped_read[0:5] 

#number of mapped bases
num_bases=mappedBases(mapped_read)
print "number of mapped bases",num_bases



############################################
#Group the mapped_reads                    #
############################################
###get the range of numbers that comprimise the mapping quality 
nums=[]
for read in mapped_read:
    nums.append(read[4])
nums=set(nums)
print "get a feel for the range of mapping qualities in this sam file", sorted(nums)
###Get the probability MAPQ=-10*log10*Pr{mapping position is wrong} 
for num in nums:
    score=inv_logit(int(num))
    print "MAPQ and probability"
    print num,score
    
group1,group2,group3=subgroups(mapped_read)

#dinuc frequency of the mapped reads

nuc=dinuc_freq(group1)
print "dinucleotide frequency of mapped reads(p<1e-3)",nuc 

#get the percentage of reads aligned need to know number of entries in fastq file 
percent=PercentReadsAligned(group1,group2,group3,reads)
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

######read transcripts processed



t1=sam_reader("/Users/samirlal/Desktop/sam/t1.sam")
t10=sam_reader("/Users/samirlal/Desktop/sam/t10.sam")
t1_2=sam_reader("/Users/samirlal/Desktop/sam/t1_2.sam")
t10_2=sam_reader("/Users/samirlal/Desktop/sam/t10_2.sam")

##get number of mapped reads printed to screen 
mapped_read=mapped_reads(t1,True)
mapped_read=mapped_reads(t10,True)
mapped_read=mapped_reads(t1_2,True)
mapped_read=mapped_reads(t10_2,True)

raw_data=raw_count_reader("/Users/samirlal/Desktop/sam/raw_counts.txt")

### Perform the normalisation methods 

rpkm1=get_RPKM(raw_data,118898,121634,136286,135102)

# write RPKM to output 
write_RPKM_data(rpkm1,"RPKM_counts.txt")

#Visualize variability among replicates using RPKM
plotreprpkm(rpkm1,"t1")
plotreprpkm(rpkm1,"t10")
plotMAreprpkm(rpkm1,"t1")
plotMAreprpkm(rpkm1,"t10")
#######################################

####Get CV 
meth1= get_cv(rpkm1,"t1")
orig=get_cv(raw_data,"t1")


####Visualise the variation (can you see how we have reduced variation possibly due to length biases and coverage biases) 

get_boxplots(meth1,orig)
plotavg_cv(meth1,orig)


#####DEexpression statistical test
plotMA(rpkm1)#Visualise MA plot

"Provide the raw count file for DE" 
"http://www.ijbcb.org/DEB/php/onlinetool.php"
      



