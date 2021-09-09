#! /usr/bin/python
#Note: this code uses python 2.7

import sys, string
from Bio import pairwise2

import numpy as np
from matplotlib import pyplot as plt

import cPickle

Dir1 = 'JPRED_random_frag_final/'
Dir2 = 'JPRED_Random_pdb_final/'

def convert(S):

    if S== ',':
        return ''
    elif S == '-':
        return 'C'
    else:
        return S

def populate_hist(hist,A1):

    for i in A1:
        hist[int(i/0.05)] += 1

def dalphabeta(p1,p2):

    return len([x for x in xrange(len(p1)) if (p1[x] == 'H' and p2[x] == 'E' or\
                                           p1[x] == 'E' and p2[x] == 'H')])

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print 'Usage: ./find_frag_abdiscrepancies_final.py JPRED_Random_frag_final/fastas.txt'
        print 'Inputs: Two directories (from Github): JPRED_random_frag_final/ and JPred_Random_pdb_final/'
        print 'Other inputs necessary for code to run: fs_xy.pik generated from compare_preds_newSS_final.py'
        print 'Outputs: Figure 3 from Mishra, et. al.'
        sys.exit()

    f = open(sys.argv[1]).read().splitlines()

    discrepancies = []

    percentSS = []
    total_disc = []
    rlengths = []
    ssamt = []

    fp = 0

    for i in f:

        d = 0

        n1 = string.split(i,'.')[0]
        n2 = string.split(n1,'_')[0]

        seq1 = open(Dir1+n1+'.fasta').read().splitlines()[1]
        pred1= open(Dir1+n1+'.jnet').read().splitlines()[1][9:]
        pred1= string.join([convert(x) for x in pred1],'')

        seq2 = open(Dir2+n2+'.fasta').read().splitlines()[1]
        pred2 = open(Dir2+n2+'.jnet').read().splitlines()[1][9:]
        pred2= string.join([convert(x) for x in pred2],'')


        a = pairwise2.align.localxs(seq1,seq2,-1,-0.5)

        if not a:
            continue

        d = dalphabeta(pred1,pred2[a[0][-2]:a[0][-1]])

        norm = min(len([x for x in pred1 if x in ['H','E']]),
                        len([x for x in pred2[a[0][-2]:a[0][-1]] if \
                             x in ['H','E']]))

        #norm = len(pred1)

        if norm != 0:
            
            percentSS.append(float(norm)/len(seq1))
            total_disc.append(float(d)/len(pred1))
            rlengths.append(len(seq1))
            ssamt.append(norm)
            

        if norm == 0:
            continue
        else:
            dNorm = float(d)/len(pred1)

        if percentSS[-1] >0.35:

            #Print false positives
            #if total_disc[-1] >=0.2:
                #print n1, len(seq1), total_disc[-1]
                #print pred1
                #print pred2[a[0][-2]:a[0][-1]]

            fp += 1

            #print n1

        discrepancies.append(dNorm)

    hist1 = [0.0]*21

    populate_hist(hist1,discrepancies)

    fs_xy = open('fs_xy.pik')
    fs_frag_ssratios = cPickle.load(fs_xy)
    fs_frag_discrepancies = cPickle.load(fs_xy)

    plt.hold(False)
    plt.plot(percentSS,total_disc,'bo',markersize=12)
    plt.hold(True)
    plt.plot(fs_frag_ssratios,fs_frag_discrepancies,'r^',markersize=12)
    plt.legend(['Random fragments','FSRs'],fontsize=12,numpoints=1)
    plt.plot([0,0.8],[0.2,0.2],'--',color=(0.75,0.75,0.75))
    #plt.plot([0.35,0.35],[0,1],'--',color=(0.75,0.75,0.75))
    plt.ylabel('Fraction helix <-> strand discrepancies')
    plt.xlabel('Fraction of residues predicted to fold into either '+ r'$\alpha$'+'-helix or '+r'$\beta$'+'-sheet')
    plt.xlim([0.35,0.8])
    plt.ylim([-0.02,0.8])

    plt.savefig('fig3.eps',dpi=200)

    print 'False positives: ', fp

