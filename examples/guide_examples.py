""" Below is code that will be run if the module is "run",
    and not just "loaded".
"""
if __name__=='__main__':
    x = Sequence('ACTGA', DNA_Alphabet, 'x')
    print "Sequence", x, "is constructed from the symbols", x.alphabet.symbols
    print "( There are", x.count('A'), "occurrences of the symbol 'A' in", x.sequence, ")"
    y = Sequence('TACGA', DNA_Alphabet, 'y')
    print "Sequence", y, "is constructed from the symbols", y.alphabet.symbols
    print
    print "( The sub-sequence 'CG' starts at index", y.find('CG'), "of", y.sequence, ")"
    print
    sm = SubstMatrix(DNA_Alphabet)
    for a in DNA_Alphabet:
        for b in DNA_Alphabet:
            if a==b: sm.set(a, b, +2) # match
            else: sm.set(a, b, -1) # mismatch
    print "Below is a substitution matrix for the alphabet", DNA_Alphabet.symbols
    print sm
    print
    aln = align(x, y, sm, -2)
    print "Below is the alignment between x and y"
    print aln
