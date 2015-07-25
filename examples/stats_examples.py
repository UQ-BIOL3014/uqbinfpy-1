
Fisher's Exact Test:
    >>> p = getFETpval(838, 159, 78, 29, False) # the one-tailed p-value (p-value=0.998 : 0.004)
    >>> print p
    >>> p2 = getFET2tail(838, 159, 78, 29, False) # the two-tailed p-value (0.006)
    >>> print p2

Wilcoxon Ranksum: (should result in a p-value around 0.0455)
    >>> p = getRSpval([4.6,5.1,5.8,6.5,4.7,5.2,6.1,7.2,4.9,5.5,6.5], [5.2,5.6,6.8,8.1,5.3,6.2,7.7,5.4,6.3,8.0])
    >>> print p

The Z-score: (should result in a z-score of around -1.78 [number of SDs from the mean] and a p-value of 0.037)
    >>> z = getZScore([12.1, 11.2, 12.3, 11.8, 11.2, 12.3, 11.1, 13.2, 12.3, 11.6, 10.8], 10.6); print z
    >>> p = f(z); print p

if __name__=='__main__':
    "A test of Ranksum: should result in a p-value around 0.0455"
    p=getRSpval([4.6,5.1,5.8,6.5,4.7,5.2,6.1,7.2,4.9,5.5,6.5], [5.2,5.6,6.8,8.1,5.3,6.2,7.7,5.4,6.3,8.0])
    print p
    "A test of Fisher's Exact test: should result in a p-value around 0.0040"
    p=getFETpval(838, 159, 78, 29, False)
    print p
    # Prob=0.002
    # p-value=0.998 : 0.004
    # two-tail=0.006
    "A test of Z-score: should result in a z-score of around -1.78 (number of SDs from the mean)"
    z = getZScore([12.1, 11.2, 12.3, 11.8, 11.2, 12.3, 11.1, 13.2, 12.3, 11.6, 10.8], 10.6)
    print z
    "A test of Pearson correlation: the example should give r = 0.988"
    r = getPearson([1, 2, 3.3, 1, 2, 3.3], [0.1, 0.2, 0.33, 0.13, 0.21, 0.3])
    print r
