###########Test Methods


t1=sam_reader("/Users/mikael/workspace/binfpy/BIOL3014/prac_5/t1.sam")

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


len_stats=length_stats(group1,group2,group3)
print len_stats

#plot the length of all three subgroups
plot_length_distrib(group1,"p<1e-3")
plot_length_distrib(group2,"1e-3<=p<1e-2")
plot_length_distrib(group3,"1e-2<=p<1")

#plot nucleotide composition along the mapped read
data=plot_base_composition(group1,'A')
data=plot_base_composition(group1,'T')
data=plot_base_composition(group1,'C')
data=plot_base_composition(group1,'G')

########read transcripts processed



t1=sam_reader("/Users/mikael/workspace/binfpy/BIOL3014/prac_5/t1.sam")
t10=sam_reader("/Users/mikael/workspace/binfpy/BIOL3014/prac_5/t10.sam")
t1_2=sam_reader("/Users/mikael/workspace/binfpy/BIOL3014/prac_5/t1_2.sam")
t10_2=sam_reader("/Users/mikael/workspace/binfpy/BIOL3014/prac_5/t10_2.sam")

##get number of mapped reads printed to screen
mapped_read=mapped_reads(t1,True)
mapped_read=mapped_reads(t10,True)
mapped_read=mapped_reads(t1_2,True)
mapped_read=mapped_reads(t10_2,True)

raw_data=raw_count_reader("/Users/mikael/workspace/binfpy/BIOL3014/prac_5/raw_counts.txt")

### Perform the normalisation methods

rpkm1=get_RPKM(raw_data,118898,121634,136286,135102)

# write RPKM to output
write_RPKM_data(rpkm1,"/Users/mikael/workspace/binfpy/BIOL3014/prac_5/RPKM_counts.txt")

#Visualize variability among replicates using RPKM
plotreprpkm(rpkm1,"t1")
plotreprpkm(rpkm1,"t10")
plotMAreprpkm(rpkm1,"t1")
plotMAreprpkm(rpkm1,"t10")
#######################################

####Get CV
meth1= get_cv(rpkm1,"t1")
orig=get_cv(raw_data,"t1")
meth2= get_cv(rpkm1,"t10")
orig2=get_cv(raw_data,"t10")

####Visualise the variation (can you see how we have reduced variation possibly due to length biases and coverage biases)

get_boxplots(meth1,orig)
plotavg_cv(meth1,orig)
get_boxplots(meth2,orig2)
plotavg_cv(meth2,orig2)


####Now try to plot MA using the FDR adjusted p-value using BH
plotMA(rpkm1)#Visualise MA plot
result_ttest=Welcht(rpkm1)
plotMA_pval(result_ttest,0.01)#plot those with corrected p-value less than 0.005


####Get diff expressed genes


sig_gene1=[]
sig_gene2=[]
for i in range(0,len(rpkm1.values())):
    fc=np.log2(float(rpkm1.values()[i][2]+1)/(rpkm1.values()[i][0]+1))
    if fc<-1.5 or fc>1.5:
        sig_gene1.append(rpkm1.keys()[i])
    else:
        pass
for i in range(0,len(rpkm1.values())):
    fc2=np.log2(float(rpkm1.values()[i][3]+1)/(rpkm1.values()[i][1]+1))
    if fc2<-1.5 or fc2>1.5:
        sig_gene2.append(rpkm1.keys()[i])
    else:
        pass

common=list(set(sig_gene1).intersection(set(sig_gene2)))

diff_express_t1={}
diff_express_t1_both={}
diff_express_t10={}
print"Genes significant by Welch t-test p<0.01"
for i in range(0,len(result_ttest)):
    if result_ttest.values()[i][4]<0.01:
        if result_ttest.keys()[i] in common:
            print result_ttest.keys()[i]
            diff_express_t1[result_ttest.keys()[i]]=result_ttest.values()[i][0] #take the first replicate
            diff_express_t10[result_ttest.keys()[i]]=result_ttest.values()[i][2]#take first replicate


t1_diff=np.array(diff_express_t1.values())
t10_diff=np.array(diff_express_t10.values())



######################check plots in current directory   #coexpression through distance based methods
#cluster the data
dend_t1=cluster_data(t1_diff,diff_express_t1.keys(),"t1")
dendt10=cluster_data(t10_diff,diff_express_t10.keys(),"t10")

#produce heatmap
heatmap_cluster(t1_diff,'t1')
heatmap_cluster(t10_diff,'t10')
