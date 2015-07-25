if __name__=='__main__': # examples to run unless this module is merely "imported"

    jp = Joint([DNA_Alphabet, DNA_Alphabet])
    jp.observe('AC', 2)
    jp.observe('A*', 3)

    pseudo = {'A':1,'C':2,'G':5,'T':3}
    d = Distrib(DNA_Alphabet, pseudo)
    d.observe('A', 3)

    ijp = IndepJoint([DNA_Alphabet, DNA_Alphabet], pseudo = 1.0)
    ijp.observe('AC')

