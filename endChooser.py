# This script chooses the last two nucleotides of allele-specific primer

import argparse
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
    print(seqRef)
    print(seqAlt)
    maxVal=0
    primerEnds=[]
    for key,item in mama[seqRef[pos-1:pos+1]][seqAlt[pos-1:pos+1]].items():
        if item>maxVal:
            primerEnds=[[key,1]]
            maxVal=item
        elif item==maxVal:
            primerEnds.append([key,1])
    print(seqRef[pos-1:pos+1],seqAlt[pos-1:pos+1])
    for key,item in mama[revComplement(seqRef[pos:pos+2])][revComplement(seqAlt[pos:pos+2])].items():
        if item>maxVal:
            primerEnds=[[key,-1]]
            maxVal=item
        elif item==maxVal:
            primerEnds.append([key,-1])
    print(revComplement(seqRef[pos:pos+2]),revComplement(seqAlt[pos:pos+2]))
    bestPrimers=[]
    for pe in primerEnds:
        if pe[1]>0:
            bestPrimers.append([seqRef[pos-30:pos-1]+pe[0],1])
        else:
            bestPrimers.append([revComplement(seqRef[pos+2:pos+30])+pe[0],-1])
    for bp in bestPrimers:
        if bp[1]>0:
            rFile.write(seqNames[-1]+'\t+\t'+bp[0]+'\t'+str(maxVal)+'\n')
        else:
            rFile.write(seqNames[-1]+'\t-\t'+bp[0]+'\t'+str(maxVal)+'\n')

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

par=argparse.ArgumentParser(description='This script chooses the last two nucleotides of allele-specific primer')
par.add_argument('--fasta-file','-ff',dest='fastaFile',type=str,help='multi-fasta-file with highlighted variable position in the format of [A/T], where A is a reference allele and T is an alternative',required=True)
args=par.parse_args()

fastaFile=args.fastaFile
try:
    f=open(fastaFile)
except FileNotFoundError:
    print('#######\nERROR! File not found:',fastaFile,'\n#######')
    exit(0)

rFile=open(fastaFile+'.primers.xls','w')
rFile.write('Sequence_Name\tStrand\tBest_Primer\tMaximal_Discrimination_Value\n')
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