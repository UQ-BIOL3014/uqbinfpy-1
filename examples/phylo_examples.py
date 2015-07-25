
""" ----------------------------------------------------------------------------------------
    Example main functions for
    (0) reading and searching tree data
    (1) running maximum parsimony to recover ancient states
    (2) running UPGMA to generate tree from alignment using evolutionary model
    ----------------------------------------------------------------------------------------"""

if __name__ == '__main__0':
    tree = parseNewick('((A:0.6,((B:3.3,(C:1.0,D:2.5)cd:1.8)bcd:5,((E:3.9,F:4.5)ef:2.5,G:0.3)efg:7)X:3.2)Y:0.5,H:1.1)I:0.2')
    print tree
    node = tree.findLabel('B')
    print node.label, node.dist, node.left, node.right
    desc = tree.getDescendantsOf(node, transitive=False)
    print 'Direct descendants of', node.label, 'are'
    for d in desc:
        print '\t', d.label
    desc = tree.getDescendantsOf(node, transitive=True)
    print 'All descendants of', node.label , 'are'
    for d in desc:
        print '\t', d.label
    a = tree.getAncestorsOf(node, transitive=False)
    print 'Direct ancestor of', node.label, 'is', a.label
    anc = tree.getAncestorsOf(node, transitive=True)
    print 'All ancestors of', node.label, 'are'
    for a in anc:
        print '\t', a.label
    from sequence import *
    tree = readNewick('/Users/mikael/workspace/binf/data/cyp3.newick')
    aln = readClustalFile('/Users/mikael/workspace/binf/data/cyp3.aln', Protein_Alphabet)
    for seq in aln:
        node = tree.findLabel(seq.name)
        a = tree.getAncestorsOf(node)
        print seq.name, 'has-distance', node.dist, 'from ancestor', a.label

if __name__ == '__main__1':
    from sequence import *
    tree = readNewick('/Users/mikael/workspace/binf/data/cyp1a1.tree')
    aln = readClustalFile('/Users/mikael/workspace/binf/data/cyp1a1.aln', Protein_Alphabet)
    tree.putAlignment(aln)
    tree.parsimony()
    print tree.root
    print tree.strSequences(10, 15)

if __name__ == '__main__2':
    from sequence import *
    tree = readNewick('/Users/mikael/Desktop/boost_tree.nwk')
    aln = readClustalFile('/Users/mikael/Desktop/interaction.aln', Bool_Alphabet)
    tree.putAlignment(aln)
    tree.parsimony()
    print tree.root
    print tree.strSequences()

if __name__ == '__main__1':
    from sequence import *
    aln = readClustalFile('/Users/mikael/workspace/binfpy/BINF6000/ws2/MalS.clustal', Protein_Alphabet)
    tree = runUPGMA(aln, 'fractional')
    print tree.root

if __name__ == '__main__':
    from sequence import *
    a = Sequence('ACG',name='a')
    b = Sequence('AAA',name='b')
    c = Sequence('TCC',name='c')
    d = Sequence('TTC',name='d')
    e = Sequence('AAA',name='e')
    aln = Alignment([a,b,c,d,e])
    D = aln.calcDistances('fractional')
    tree = runUPGMA(aln,'fractional')
    tree.putAlignment(aln)
    print tree
    tree.parsimony()
    print tree.root
    print tree.strSequences(0, 3)
