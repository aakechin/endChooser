# This script chooses the last two nucleotides of allele-specific primer

import argparse
import glob
from Bio.Seq import Seq
import os
import primer3
from operator import itemgetter
from Bio.SeqUtils import MeltingTemp as mt

thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

def revComplement(nuc):
    return(str(Seq(nuc).reverse_complement()))

def constructPrimers(mutPrimerSeq,wtPrimerSeq,seq,ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg):
    # First of all choose the left primers. It's got already defined 3'-end
    primers=[]
    primerProps={}
    primerScores={}
    for i in range(primerLen-primerLenDev,primerLen+primerLenDev+1):
        mutPrimer=mutPrimerSeq[0][-i:]
        tm1=mt.Tm_Wallace(mutPrimer)
        homodG1=primer3.calcHomodimer(mutPrimer).dg
        hairpindG1=primer3.calcHairpin(mutPrimer).dg
        for z in range(primerLen-primerLenDev,primerLen+primerLenDev+1):
            wtPrimer=wtPrimerSeq[-z:]
            tm2=mt.Tm_Wallace(wtPrimer)
            homodG2=primer3.calcHomodimer(wtPrimer).dg
            hairpindG2=primer3.calcHairpin(wtPrimer).dg
            heterodG0=primer3.calcHeterodimer(mutPrimer,wtPrimer).dg
            for j in range(ampliconLen-ampliconLenDev,ampliconLen+ampliconLenDev+1):
                rPrimerStart=seq.index(wtPrimer)+j
                if rPrimerStart>len(seq):
                    print('ERROR! The length of input sequence is not enough for construction of R-primer!')
                    print('Input sequence length:',len(seq))
                    print('Position of substitution:',seq.index(wtPrimer)-len(wtPrimer))
                    print('Maximal amplicon length:',ampliconLen+ampliconLenDev)
                    exit(0)
                for k in range(primerLen-primerLenDev,primerLen+primerLenDev+1):
                    rPrimer=revComplement(seq[rPrimerStart-k:rPrimerStart+1])
                    tm3=mt.Tm_Wallace(rPrimer)
                    homodG3=primer3.calcHomodimer(rPrimer).dg
                    heterodG1=primer3.calcHeterodimer(rPrimer,mutPrimer).dg
                    heterodG2=primer3.calcHeterodimer(rPrimer,wtPrimer).dg
                    hairpindG3=primer3.calcHairpin(rPrimer).dg
                    primers.append([mutPrimer,wtPrimer,rPrimer])
                    primerProps['_'.join(primers[-1])]=[len(mutPrimer),len(wtPrimer),len(rPrimer),tm1,tm2,tm3,homodG1,homodG2,homodG3,hairpindG1,hairpindG2,hairpindG3,heterodG0,heterodG1,heterodG2]
                    primerScores['_'.join(primers[-1])]=20*(abs(tm1-meltTemp)+abs(tm2-meltTemp)+abs(tm3-meltTemp)+abs(tm1-tm2)+abs(tm2-tm3)+abs(tm1-tm3))+5*(min(homodG1,dimerdg)+min(homodG2,dimerdg)+min(homodG3,dimerdg)+min(heterodG0,dimerdg)+min(heterodG1,dimerdg)+min(heterodG2,0))/(6*dimerdg)+(min(hairpindG1,dimerdg)+min(hairpindG2,dimerdg)+min(hairpindG3,dimerdg))/(3*dimerdg)
    theBest=None
    bestMatch=None
    for key,value in sorted(primerScores.items(),key=itemgetter(1),reverse=False):
        if theBest is None:
            theBest=[key,value,primerProps[key]]
        if ((primerProps[key][3]<=meltTemp+meltTempDev and primerProps[key][4]<=meltTemp+meltTempDev and primerProps[key][5]<=meltTemp+meltTempDev
            and primerProps[key][3]>=meltTemp-meltTempDev and primerProps[key][4]>=meltTemp-meltTempDev and primerProps[key][5]>=meltTemp-meltTempDev
            and primerProps[key][6]>=dimerdg and primerProps[key][7]>=dimerdg and primerProps[key][8]>=dimerdg
            and primerProps[key][9]>=dimerdg and primerProps[key][10]>=dimerdg and primerProps[key][11]>=dimerdg
            and primerProps[key][12]>=dimerdg and primerProps[key][13]>=dimerdg and primerProps[key][14]>=dimerdg) and bestMatch is None):
            primerProps[key]=map(str,primerProps[key])
            bestMatch=[key,value,primerProps[key]]
    theBest[2]=map(str,theBest[2])
    return(bestMatch,theBest)

def chooseBestPrimers(seq,seqNames,rFile,ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg):
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
    wtPrimers=[]
    bestMatchPrimers=[]
    for pe in primerEnds:
        if pe[1]>0:
            bestPrimers.append([seqRef[max(pos-primerLen-primerLenDev,0):pos-1]+pe[0],1])
            wtPrimers.append(seqRef[max(pos-primerLen-primerLenDev,0):pos+1])
            bestMatch,theBest=constructPrimers(bestPrimers[-1],wtPrimers[-1],seqRef,ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg)
            bestMatchPrimers.append([bestMatch,theBest])
            plus40.append(seqRef[pos+1:min(len(seqRef),pos+primerLen+primerLenDev+1)])
        else:
            bestPrimers.append([revComplement(seqRef[pos+2:min(len(seqRef),pos+primerLen+primerLenDev)])+pe[0],-1])
            wtPrimers.append(revComplement(seqRef[pos:min(len(seqRef),pos+primerLen+primerLenDev)]))
            bestMatch,theBest=constructPrimers(bestPrimers[-1],wtPrimers[-1],revComplement(seqRef),ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg)
            bestMatchPrimers.append([bestMatch,theBest])
            plus40.append(revComplement(seqRef[max(pos-primerLen-primerLenDev,0):pos]))
    for i,bp in enumerate(bestPrimers):
        if bp[1]>0:
            if bestMatchPrimers[i][0]:
                rFile.write(seqNames[-1]+'\t'+ref+'/'+alt+'\t'+'+\t'+'\t'.join(bestMatchPrimers[i][0][0].split('_'))+'\t'+str(primerVals[i])+'\t'+'\t'.join(bestMatchPrimers[i][0][2])+'\n')
            else:
                rFile.write(seqNames[-1]+'\t'+ref+'/'+alt+'\t'+'+\t'+'\t'.join(bestMatchPrimers[i][1][0].split('_'))+'\t'+str(primerVals[i])+'\t'+'\t'.join(bestMatchPrimers[i][1][2])+'\n')
        else:
            if bestMatchPrimers[i][0]:
                rFile.write(seqNames[-1]+'\t'+revComplement(ref)+'/'+revComplement(alt)+'\t'+'-\t'+'\t'.join(bestMatchPrimers[i][0][0].split('_'))+'\t'+str(primerVals[i])+'\t'+'\t'.join(bestMatchPrimers[i][0][2])+'\n')
            else:
                rFile.write(seqNames[-1]+'\t'+revComplement(ref)+'/'+revComplement(alt)+'\t'+'-\t'+'\t'.join(bestMatchPrimers[i][1][0].split('_'))+'\t'+str(primerVals[i])+'\t'+'\t'.join(bestMatchPrimers[i][1][2])+'\n')

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
par.add_argument('--amplicon-length','-al',dest='ampliconLen',type=str,help='average length of amplicon. Default: 150 bp',default=150,required=False)
par.add_argument('--amplicon-length-deviation','-ald',dest='ampliconLenDev',type=str,help='deviation of amplicon length. Default: 10 bp',default=10,required=False)
par.add_argument('--primer-length','-pl',dest='primerLen',type=str,help='optimal primer length. Default: 20 nucleotides',default=20,required=False)
par.add_argument('--primer-length-deviation','-pld',dest='primerLenDev',type=str,help='deviation of primer length. Default: 4 nucleotides',default=4,required=False)
par.add_argument('--melt-temp','-at',dest='meltTemp',type=str,help='optimal temperature of annealing. Default: 60 degrees Celsius',default=60,required=False)
par.add_argument('--melt-temp-deviation','-atd',dest='meltTempDev',type=str,help='deviation of temperature of annealing. Default: 5 degrees Celsius',default=5,required=False)
par.add_argument('--min-dimer-dg','-dimerdg',dest='dimerdg',type=str,help='minimal value of dG of dimer formation. Default: -3000 kcal/mol',default=-3000,required=False)
args=par.parse_args()

fastaFile=args.fastaFile
ampliconLen=args.ampliconLen
ampliconLenDev=args.ampliconLenDev
primerLen=args.primerLen
primerLenDev=args.primerLenDev
meltTemp=args.meltTemp
meltTempDev=args.meltTempDev
dimerdg=args.dimerdg
try:
    f=open(fastaFile)
except FileNotFoundError:
    print('#######\nERROR! File not found:',fastaFile,'\n#######')
    exit(0)

rFile=open(fastaFile+'.primers.xls','w')
rFile.write('Sequence_Name\tRef/Alt\tStrand\tBest_Primer_for_Mutant\tBest_Primer_for_WT\tBest_Reverse_Primer\tDiscrimination_Value\tLen1\tLen2\tLen3\tTm1\tTm2\tTm3\tHomodimer1\tHomodimer2\tHomodimer3\tHairpin1\tHairpin2\tHairpin3\tHeterodimer1\tHeterodimer2\tHeterodimer3\n')
seqNames=[]
for line in f:
    if '>' in line:
        if len(seqNames)>0:
            chooseBestPrimers(seq,seqNames,rFile,ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg)
        seqNames.append(line.replace('> ','').replace('>','').replace('\n','').replace('\r',''))
        seq=''
    else:
        seq+=line.replace('\n','').replace('\r','').replace(' ','').replace('\t','')
chooseBestPrimers(seq,seqNames,rFile,ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg)
rFile.close()
print('Done')
