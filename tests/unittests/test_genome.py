import sys
import os
import pytest

sys.path.insert(0, os.path.abspath('../../'))

from uqbinfpy import genome

def test_all():
    # ------------------- Example ---------------------

    ge3716 = readGEOFile('/Users/mikael/workspace/COSC2000/GDS3716.soft')

    ratio = GeneExpression('GDS3716_ratio')
    ratio.addSamples('S1_ER+/Healthy', ge3716.getRatio( 33,  0))
    ratio.addSamples('S2_ER+/Healthy', ge3716.getRatio( 34,  1))
    ratio.addSamples('S3_ER+/Healthy', ge3716.getRatio( 35,  2))
    ratio.addSamples('S4_ER+/Healthy', ge3716.getRatio( 36,  3))
    ratio.addSamples('S5_ER+/Healthy', ge3716.getRatio( 37,  4))
    ratio.addSamples('S6_ER+/Healthy', ge3716.getRatio( 38,  5))
    ratio.addSamples('S7_ER+/Healthy', ge3716.getRatio( 39,  6))
    ratio.addSamples('S8_ER+/Healthy', ge3716.getRatio( 40,  7))
    ratio.addSamples('S9_ER+/Healthy', ge3716.getRatio( 41,  8))
    ratio.addSamples('S1_ER-/Healthy', ge3716.getRatio( 24,  9))
    ratio.addSamples('S2_ER-/Healthy', ge3716.getRatio( 25, 10))
    ratio.addSamples('S3_ER-/Healthy', ge3716.getRatio( 26, 11))
    ratio.addSamples('S4_ER-/Healthy', ge3716.getRatio( 27, 12))
    ratio.addSamples('S5_ER-/Healthy', ge3716.getRatio( 28, 13))
    ratio.addSamples('S6_ER-/Healthy', ge3716.getRatio( 29, 14))
    ratio.addSamples('S7_ER-/Healthy', ge3716.getRatio( 30, 15))
    ratio.addSamples('S8_ER-/Healthy', ge3716.getRatio( 31, 16))
    ratio.addSamples('S9_ER-/Healthy', ge3716.getRatio( 32, 17))
    ratio.writeGEOFile('/Users/mikael/workspace/COSC2000/GDS3716_ratios.soft')
    print ge3716.getHeaders()


    z = ratio.getZScore(0) # NOT recommended! Ratios are NOT normally distributed! Use log-ratios instead.

    ge38 = readGEOFile('/Users/mikael/workspace/COSC2000/GDS38.soft', id_column = 1)
    cln2_profile = ge38.getGenes('CLN2')
    pcorr = ge38.getPearson('CLN2')
    gp = GeneExpression('Ex3', 'PC_CLN2', pcorr)
    sorted = gp.sort('PC_CLN2', True)
    print sorted[0], ge38.getGenes(sorted[0])
    print sorted[1], ge38.getGenes(sorted[1])

