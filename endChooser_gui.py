# This script chooses the last two nucleotides of allele-specific primer

import easygui
import glob
from Bio.Seq import Seq
import os

thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

def revComplement(nuc):
    return(str(Seq(nuc).reverse_complement()))

def chooseBestPrimers(seq,seqNames,rFile):
    ref=seq[seq.find('[')+1:seq.find('/')]
    alt=seq[seq.find('/')+1:seq.find(']')]
    pos=seq.find('[')
    seqRef=seq.replace(seq[seq.find('['):seq.find(']')+1],ref)
    seqAlt=seq.replace(seq[seq.find('['):seq.find(']')+1],alt)
    maxVal=0
    primerVals=[]
    primerEnds=[]
    for key,item in mama[seqRef[pos-1:pos+1]][seqAlt[pos-1:pos+1]].items():
        primerEnds.append([key,1])
        primerVals.append(item)
    for key,item in mama[revComplement(seqRef[pos:pos+2])][revComplement(seqAlt[pos:pos+2])].items():
        primerEnds.append([key,-1])
        primerVals.append(item)
    bestPrimers=[]
    plus40=[]
    for pe in primerEnds:
        if pe[1]>0:
            bestPrimers.append([seqRef[max(pos-40,0):pos-1]+pe[0],1])
            plus40.append(seqRef[pos+1:min(len(seqRef),pos+41)])
        else:
            bestPrimers.append([revComplement(seqRef[pos+2:min(len(seqRef),pos+40)])+pe[0],-1])
            plus40.append(revComplement(seqRef[max(pos-40,0):pos]))
    for i,bp in enumerate(bestPrimers):
        if bp[1]>0:
            rFile.write(seqNames[-1]+'\t'+ref+'/'+alt+'\t'+'+\t'+bp[0]+'\t'+str(primerVals[i])+'\t'+plus40[i]+'\n')
        else:
            rFile.write(seqNames[-1]+'\t'+revComplement(ref)+'/'+revComplement(alt)+'\t'+'-\t'+bp[0]+'\t'+str(primerVals[i])+'\t'+plus40[i]+'\n')

mama={}
file=open(thisDir+'mamaPrimers.txt')
for string in file:
    if string=='' or string=='\n': continue
    cols=string.replace('\n','').replace('\r','').split('\t')
    if cols[0] not in mama.keys():
        mama[cols[0]]={}
    if cols[1] not in mama[cols[0]].keys():
        mama[cols[0]][cols[1]]={}
    alts=cols[2].split('/')
    for alt in alts:
        mama[cols[0]][cols[1]][alt]=int(cols[3])            

fastaFile=easygui.fileopenbox()

try:
    f=open(fastaFile)
except FileNotFoundError:
    print('#######\nERROR! File not found:',fastaFile,'\n#######')
    exit(0)

rFile=open(fastaFile+'.primers.xls','w')
rFile.write('Sequence_Name\tRef/Alt\tStrand\tBest_Primer\tMaximal_Discrimination_Value\t+40bp_after_primer\n')
seqNames=[]
for line in f:
    if '>' in line:
        if len(seqNames)>0:
            chooseBestPrimers(seq,seqNames,rFile)
        seqNames.append(line.replace('> ','').replace('>','').replace('\n','').replace('\r',''))
        seq=''
    else:
        seq+=line.replace('\n','').replace('\r','')
chooseBestPrimers(seq,seqNames,rFile)
rFile.close()
