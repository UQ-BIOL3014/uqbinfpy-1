Examples:

>>> import seqdata as sqd
>>> hg19 = sqd.TwoBitFile('/Users/mikael/simhome/hg19.2bit')
>>> mybed1 = sqd.BedFile('/Users/mikael/simhome/mcf7/h3k4me3.bed', format = 'Peaks')
>>> print mybed1[0]
    ('chr1', 540640, 540790)
>>> for e in mybed1[0:5]:
...    print e.getloc(hg19)

    ggcgttttcctgtaaagttgggcacacgcttcccacatgactcagcaattgcacttctgggtatgtacccgagagaaacaaaagcttatgttcacacaaaaacctacaacgcaaatgcacaaacagctctatccaacaaccctggaagca
    ATATAGTAAAACCCAGCCCATGGCCCCTAACAGGGGCCCTCTCAGCCCTCCTAATGACCTCCGGCCTAGCCATGTGATTTCACTTCCACTCCACAACCCTCCTCATACTAGGCCTACTAACCAACACACTAACCATATACCAATGATGGC
    ...
>>> for key in hg19:
        print key

    chr19_gl000208_random
    chr8_gl000197_random
    chr6_apd_hap1
    chr13
    ...
>>> hg19['chrX']
    Out[1]: <seqdata.TwoBitSequence at 0x10cde7d90>
>>> len(hg19['chrX'])
    Out[1]: 155270560
>>> hg19['chrX'][1000000:1000060]
    Out[1]: 'AAAcagctacttggaaggctgaagcaggaggattgtttgagtctaggagtttgaggctgc'

