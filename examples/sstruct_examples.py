""" -------------------------------------------
Below is test code
------------------------------------------- """

# Read some protein sequence data
#prot = sequence.readFastaFile('/Users/mikael/workspace/binf/data/prot2.fa', prot_alpha)
# read the secondary structure data for the proteins above (indices should agree)
#sstr = sequence.readFastaFile('/Users/mikael/workspace/binf/data/sstr3.fa', sstr_alpha)

prot = [sequence.Sequence('PNKRKGFSEGLWEIENNPTVKASGY', symbol.Protein_Alphabet, '2NLU_r76')]
sstr = [sequence.Sequence('CCCCHHHHHHHHHHHCCCCCCCCCC', symbol.DSSP3_Alphabet, '2NLU_s76')]

tp = 0 # number of true positives (correctly identified calls)
tn = 0 # number of true negatives (correctly missed no-calls)
fp = 0 # number of false positives (incorrectly identified no-calls)
fn = 0 # number of false negatives (incorrectly missed calls)

for index in range(len(prot)):

    myprot = prot[index]
    mysstr = sstr[index]
    myalpha = [sym == 'H' for sym in sstr[index]]
    mybeta = [sym == 'E' for sym in sstr[index]]

    """
     1. Assign all of the residues in the peptide the appropriate set of parameters.
    """
    alpha = getScores(myprot, 0)
    beta  = getScores(myprot, 1)
    turn  = getScores(myprot, 2)

    """
     2. Scan through the peptide and identify regions where 4 out of 6 contiguous residues have P(a-helix) > 100.
        That region is declared an alpha-helix.
        Extend the helix in both directions until a set of four contiguous residues that have an average P(a-helix) < 100 is reached.
        That is declared the end of the helix.
        If the segment defined by this procedure is longer than 5 residues and the average P(a-helix) > P(b-sheet) for that segment,
        the segment can be assigned as a helix.

     3. Repeat this procedure to locate all of the helical regions in the sequence.
    """
    calls_a1 = markCountAbove(alpha, width = 6, call_cnt = 4)
    calls_a2 = extendDownstream(alpha, calls_a1, width = 4)
    calls_a3 = extendUpstream(alpha, calls_a2, width = 4)

    """
     4. Scan through the peptide and identify a region where 3 out of 5 of the residues have a value of P(b-sheet) > 100.
        That region is declared as a beta-sheet.
        Extend the sheet in both directions until a set of four contiguous residues that have an average P(b-sheet) < 100 is reached.
        That is declared the end of the beta-sheet.
        Any segment of the region located by this procedure is assigned as a beta-sheet
        if the average P(b-sheet) > 105 and the average P(b-sheet) > P(a-helix) for that region.
    """
    calls_b1 = markCountAbove(beta, width = 5, call_cnt = 3)
    calls_b2 = extendDownstream(beta, calls_b1, width = 4)
    calls_b3 = extendUpstream(beta, calls_b2, width = 4)

    """
     5. Any region containing overlapping alpha-helical and beta-sheet assignments are taken to be helical
        if the average P(a-helix) > P(b-sheet) for that region.
        It is a beta sheet if the average P(b-sheet) > P(a-helix) for that region.
    """
    avg_a = calcRegionAverage(alpha, calls_a3)
    avg_b = calcRegionAverage(beta, calls_b3)
    diff_a = [avg_a[i] - avg_b[i] for i in range(len(avg_a))]
    diff_b = [avg_b[i] - avg_a[i] for i in range(len(avg_a))]
    calls_a4 = checkSupport(calls_a3, diff_a)
    calls_b4 = checkSupport(calls_b3, diff_b)

    i = 0
    for call in myalpha:
        if call:
            if calls_a4[i]: tp += 1
            else: fn += 1
        else:
            if calls_a4[i]: fp += 1
            else: tn += 1
        i += 1

print "TP = %d" % tp
print "TN = %d" % tn
print "FP = %d" % fp
print "FN = %d" % fn
print "Accuracy = %d%%" % ((tp + tn) * 100 / (tp + tn + fp + fn))

print myprot
print mysstr
print 'P(Helix):', makesstr(calls_a4, 'H')
print 'P(Sheet):', makesstr(calls_b4, 'E')
