if __name__ == '__main__0':
    writeBitFile('/Users/mikael/simhome/share/gene_association.tair',
                 '/Users/mikael/simhome/share/gene_ontology_ext.obo',
                 '/Users/mikael/simhome/share/tair.bit')
    """
    Started at Tue Jul 17 09:24:25 2012
    Read annotations for 14973326 genes
    Read 35980 GO definitions
    Wrote header 14973326    35980    6    7    19, now at @302
    Wrote GO annotations, now at @773175002
    Wrote 35980 GO definitions, now at @775623250
    Completed at Tue Jul 17 10:59:05 2012
    """
    bgo = BinGO('/Users/mikael/simhome/TNR/GO/tair.bit')
    #bgo = BinGO('/Users/mikael/simhome/gene_association.bit', taxa = 39947)
    print "Done loading index with %d genes annotated" % len(bgo.annot_index)

if __name__ == '__main__1':
    print os.getcwd()
    #writeBitFile('/Users/mikael/simhome/gene_association.goa_uniprot',
    #             '/Users/mikael/simhome/gene_ontology_ext.obo',
    #             '/Users/mikael/simhome/gene_association_mammal.bit',
    ##             [9606,10090])
    """
    Started at Tue Jul 17 09:24:25 2012
    Read annotations for 14973326 genes
    Read 35980 GO definitions
    Wrote header 14973326    35980    6    7    19, now at @302
    Wrote GO annotations, now at @773175002
    Wrote 35980 GO definitions, now at @775623250
    Completed at Tue Jul 17 10:59:05 2012
    """
    bgo = BinGO('/Users/mikael/simhome/share/goa.bit')
    #bgo = BinGO('/Users/mikael/simhome/gene_association.bit', taxa = 39947)
    print "Done loading index with %d genes annotated" % len(bgo.annot_index)
    #pos = [id.strip() for id in open('/Users/mikael/simhome/nls/identifiers_streptophyta_pos.txt')]
    #neg = [id.strip() for id in open('/Users/mikael/simhome/nls/identifiers_streptophyta_neg.txt')]
    #rows = bgo.getGOReport(pos, background = neg, taxa = 39947, include_more_general = True)
    f_bg = open('/Users/mikael/simhome/homoaa/uniprotID/S288C-sgd.proteinID.list')
    f_fg = open('/Users/mikael/simhome/homoaa/uniprotID/S288C-sgd.homo')
    bg = [s.strip() for s in f_bg]
    print 'Background has %d proteins' % len(bg)
    fg = []
    for tract in f_fg:
        field = tract.split('\t')
        id = field[1]
        aa_cnt = int(field[5])
        tn_cnt = int(field[6])
        if aa_cnt >= 7:
            fg.append(id)
    print 'Foreground has %d proteins' % len(fg)
    f_bg.close()
    f_fg.close()
    rows = bgo.getGOReport(fg, bg, include_more_general = True)
    for row in rows[0:100]:
        if len(row) > 4:
            print "%s\t%4.2E\t%3d\t%6d\t%s (%s)" % (row[0], row[1], row[2], row[3], row[4].strip(), row[5])
        else:
            print "%s\t%3d\t%s (%s)" % (row[0], row[1], row[2].strip(), row[3])

if __name__ == '__main__':
    go = GO('/Users/mikael/simhome/share/gene_association.goa_ref_human', '/Users/mikael/simhome/share/go-basic.obo')
    myproteins = ['Q9H9L7', 'P41223', 'Q13352', 'Q9UFW8', 'P23528', 'P21291', 'P50461', 'P60981', 'P07992', 'P51858', 'Q0VD86', 'P61244', 'O60682', 'Q01658', 'Q9HAN9', 'P06748', 'P52945', 'Q6MZT1', 'P61571', 'Q8NHV9', 'Q96EU6', 'Q14493', 'Q9NS25', 'Q8IZU3', 'Q9P016', 'Q96B42', 'O60688']
    mybackground = go.getGenes4Term('GO:0005634', include_more_specific = True)
    print len(mybackground)
    report = go.getEnrichmentReport(myproteins, mybackground, threshold = 0.05, include_more_general = True)
    for row in report:
        print row
