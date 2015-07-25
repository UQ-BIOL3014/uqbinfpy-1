if __name__ == '__main__':
    aln = readClustalFile('/Users/mikael/simhome/ASR/dp16_example.aln', Protein_Alphabet)
    aln.displayConsensus(theta1 = 0.1, theta2 = 0.01)
