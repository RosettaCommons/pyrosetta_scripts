#! /usr/bin/env python
"""
Created on Fri Dec  2 14:37:17 2016

@author: hamedkhakzad
"""

## This function helps to extract each XL peptide from sequence.
## For PEPKPEP--QPEPKQPEP we need to call this function two times to give us
## both peptides. It also helps to define position of K on each peptide which
## is needed in kojak format.
def kojak_generator (sequence, pos):

    xlink = []
    pos_on_pep = 0

    if sequence[pos] == 'K':
        pos_temp1=0
        pos_temp2=0
        
        for cnt1 in range(pos-1, -1, -1):
            if (sequence[cnt1] == 'K') or (sequence[cnt1] == 'R'):
                pos_temp1 = cnt1 + 1
                break
            else:
                pos_temp1 = 0
        
        for cnt2 in range(pos+1, len(sequence), +1):
            if (sequence[cnt2] == 'K') or (sequence[cnt2] == 'R'):
                pos_temp2 = cnt2 + 1
                break
            else:
                pos_temp2 = len(sequence)
         
        xlink = sequence[pos_temp1:pos_temp2]
        pos_on_pep = pos - pos_temp1 + 1
    else:
        print "ERROR"
    
    return xlink , pos_on_pep