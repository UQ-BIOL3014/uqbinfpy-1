# Example 1: Find the peroxisome targeting signal
if __name__=='__main__0':
    import os
    os.chdir('/Users/mikael/workspace/binf/data/')  # set to the directory where you keep your files
    seqs = sequence.readFastaFile('pex2.fa', symbol.Protein_Alphabet)
    W = 3
    pseudo = prob.readDistrib('blosum62.distrib')
    gibbs = GibbsMotif(seqs, W)
    q = gibbs.discover(pseudo)
    p = gibbs.getBackground()

    # Let's display the results, i.e. the best matches to the found motif
    a = getAlignment(seqs, q, p)
    k = 0
    for seq in seqs:
        print "%s \t%d \t%s" % (str(seq), a[k], seq[a[k]:a[k]+W])
        k += 1

    # save the motif in two files: one for the foreground distributions and one with the background
    prob.writeDistribs(q, 'pex2q.distrib')
    p.writeDistrib('pex2p.distrib')

    # the data can be re-loaded and used the following way
    new_q = prob.readDistribs('pex2q.distrib')
    new_p = prob.readDistrib('pex2p.distrib')
    a = getAlignment(seqs, new_q, new_p) # the best alignment is identified from loaded motif
    gibbs2 = GibbsMotif(seqs, W, a) # alignment is provided to bootstrap training
    gibbs2.discover(pseudo, 20000)  # train

# Example 2: Find DNA binding signals in some chip-seq data
# This data set is large so takes time... ~10 mins on my quad-core iMac; 100000 rounds are sufficient in my tests, should reach LL~9000
if __name__=='__main__1':
    import os
    os.chdir('/Users/mikael/workspace/binf/data/')  # set to the directory where you keep your files
    seqs = sequence.readFastaFile('chipseq_2330.fa', symbol.DNA_Alphabet)
    W = 10
    gibbs = GibbsMotif(seqs, W)
    q = gibbs.discover(niter = 100000)
    p = gibbs.getBackground()
    a = getAlignment(seqs, q, p)
    k = 0
    for seq in seqs:
        print "%s \t%d \t%s" % (seq.name, a[k], seq[a[k]:a[k]+W])
        k += 1
    prob.writeDistribs(q, 'chipseq_2330_W10q.distrib')
    p.writeDistrib('chipseq_2330_W10p.distrib')
