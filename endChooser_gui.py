# This script chooses the last two nucleotides of allele-specific primer

import argparse,re
import glob
from Bio.Seq import Seq
import os
import primer3
from operator import itemgetter
from Bio.SeqUtils import MeltingTemp as mt
import subprocess as sp

thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

def revComplement(nuc):
    return(str(Seq(nuc).reverse_complement()))

def constructPrimers(mutPrimerSeq,wtPrimerSeq,seq,ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg,needMinAmpl):
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
            for j in range(ampliconLen-ampliconLenDev,ampliconLen+ampliconLenDev+1):
                rPrimerStart=seq.index(wtPrimer[:-2])+j
                if rPrimerStart>len(seq):
                    print('ERROR! The length of input sequence is not enough for construction of R-primer!')
                    print('Input sequence length:',len(seq))
                    print('Position of substitution:',seq.index(wtPrimer)-len(wtPrimer))
                    print('Maximal amplicon length:',ampliconLen+ampliconLenDev)
                    exit(0)
                for k in range(primerLen-primerLenDev,primerLen+primerLenDev+1):
                    rPrimer=revComplement(seq[rPrimerStart-k:rPrimerStart+1])
                    tm3=mt.Tm_Wallace(rPrimer)
                    try:
                        homodG3=primer3.calcHomodimer(rPrimer).dg
                    except OSError:
                        homodG3=0
                    try:
                        heterodG1=primer3.calcHeterodimer(rPrimer,mutPrimer).dg
                    except OSError:
                        heterodG1=0
                    try:
                        heterodG2=primer3.calcHeterodimer(rPrimer,wtPrimer).dg
                    except OSError:
                        heterodG2=0
                    try:
                        hairpindG3=primer3.calcHairpin(rPrimer).dg
                    except OSError:
                        hairpindG3=0
                    primers.append([mutPrimer,wtPrimer,rPrimer])
                    wtPrimerStart=seq.index(wtPrimer[:-2])
                    mutPrimerStart=seq.index(mutPrimer[:-2])
                    primerProps['_'.join(primers[-1])]=[len(mutPrimer),len(wtPrimer),len(rPrimer),j,tm1,tm2,tm3,homodG1,homodG2,homodG3,hairpindG1,hairpindG2,hairpindG3,heterodG1,heterodG2,mutPrimerStart,wtPrimerStart,rPrimerStart]
                    if needMinAmpl:
                        primerScores['_'.join(primers[-1])]=j+20*(abs(tm1-meltTemp)+abs(tm2-meltTemp)+abs(tm3-meltTemp)+abs(tm1-tm2)+abs(tm2-tm3)+abs(tm1-tm3))+5*(min(homodG1,dimerdg)+min(homodG2,dimerdg)+min(homodG3,dimerdg)+min(heterodG1,dimerdg)+min(heterodG2,dimerdg))/(5*dimerdg)+(min(hairpindG1,dimerdg)+min(hairpindG2,dimerdg)+min(hairpindG3,dimerdg))/(3*dimerdg)
                    else:
                        primerScores['_'.join(primers[-1])]=20*(abs(tm1-meltTemp)+abs(tm2-meltTemp)+abs(tm3-meltTemp)+abs(tm1-tm2)+abs(tm2-tm3)+abs(tm1-tm3))+5*(min(homodG1,dimerdg)+min(homodG2,dimerdg)+min(homodG3,dimerdg)+min(heterodG1,dimerdg)+min(heterodG2,dimerdg))/(5*dimerdg)+(min(hairpindG1,dimerdg)+min(hairpindG2,dimerdg)+min(hairpindG3,dimerdg))/(3*dimerdg)
    theBest=None
    bestMatch=None
    for key,value in sorted(primerScores.items(),key=itemgetter(1),reverse=False):
        if theBest is None:
            theBest=[key,value,primerProps[key]]
        if ((primerProps[key][4]<=meltTemp+meltTempDev and primerProps[key][5]<=meltTemp+meltTempDev and primerProps[key][6]<=meltTemp+meltTempDev
            and primerProps[key][4]>=meltTemp-meltTempDev and primerProps[key][5]>=meltTemp-meltTempDev and primerProps[key][6]>=meltTemp-meltTempDev
            and primerProps[key][7]>=dimerdg and primerProps[key][8]>=dimerdg and primerProps[key][9]>=dimerdg
            and primerProps[key][10]>=dimerdg and primerProps[key][11]>=dimerdg and primerProps[key][12]>=dimerdg
            and primerProps[key][13]>=dimerdg and primerProps[key][14]>=dimerdg) and bestMatch is None):
            primerProps[key]=map(str,primerProps[key])
            bestMatch=[key,value,primerProps[key]]
    theBest[2]=map(str,theBest[2])
    return(bestMatch,theBest)

def chooseBestPrimers(seq,seqNames,rFile,ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg,needMinAmpl):
    if '[' in seq:
        if ']' in seq:
            lBracket='['
            rBracket=']'
        else:
            print('#######\nERROR! There is incorrect designation of alternative alleles!\n'
                  'It should lool like [A/G] or (A/G)',seq[seq.find(lBracket):seq.find(rBracket)+1],'\n#######')
            exit(0)
    elif '(' in seq:
        if ')' in seq:
            lBracket='('
            rBracket=')'
        else:
            print('#######\nERROR! There is incorrect designation of alternative alleles!\n'
                  'It should lool like [A/G] or (A/G)',seq[seq.find(lBracket):seq.find(rBracket)+1],'\n#######')
            exit(0)
    else:
        print('#######\nERROR! There is incorrect designation of alternative alleles!\n'
              'It should lool like [A/G] or (A/G)',seq[seq.find(lBracket):seq.find(rBracket)+1],'\n#######')
        exit(0)
    ref=seq[seq.find(lBracket)+1:seq.find('/')]
    alt=seq[seq.find('/')+1:seq.find(rBracket)]
    if len(ref)!=1 or len(alt)!=1:
        print('#######\nERROR! There is incorrect designation of alternative alleles!\n'
              'It should lool like [A/G] or (A/G). But now it is ',seq[seq.find(lBracket):seq.find(rBracket)+1],'\n#######')
        exit(0)
    pos=seq.find(lBracket)
    seqRef=seq.replace(seq[seq.find(lBracket):seq.find(rBracket)+1],ref)
    seqAlt=seq.replace(seq[seq.find(lBracket):seq.find(rBracket)+1],alt)
    maxVal=0
    primerVals=[]
    primerEnds=[]
    # MAMA-primer on the plus strand for mutant allele
    for key,item in mama[seqRef[pos-1:pos+1]][seqAlt[pos-1:pos+1]].items():
        primerEnds.append([key,1])
        primerVals.append(item)
    # MAMA-primer on the minus strand for mutant allele
    for key,item in mama[revComplement(seqRef[pos:pos+2])][revComplement(seqAlt[pos:pos+2])].items():
        primerEnds.append([key,-1])
        primerVals.append(item)
    wtPrimerVals=[]
    wtPrimerEnds=[]
    # MAMA-primer on the plus strand for mutant allele
    for key,item in mama[seqAlt[pos-1:pos+1]][seqRef[pos-1:pos+1]].items():
        wtPrimerEnds.append([key,1])
        wtPrimerVals.append(item)
    # MAMA-primer on the minus strand for mutant allele
    for key,item in mama[revComplement(seqAlt[pos:pos+2])][revComplement(seqRef[pos:pos+2])].items():
        wtPrimerEnds.append([key,-1])
        wtPrimerVals.append(item)
    bestPrimers=[]
    plus40=[]
    wtPrimers=[]
    bestMatchPrimers=[]
    for pe,wtpe in zip(primerEnds,wtPrimerEnds):
        if pe[1]>0: # if strand is plus
            bestPrimers.append([seqRef[max(pos-primerLen-primerLenDev,0):pos-1]+pe[0],1])
            wtPrimers.append(seqRef[max(pos-primerLen-primerLenDev,0):pos-1]+wtpe[0])
            bestMatch,theBest=constructPrimers(bestPrimers[-1],wtPrimers[-1],seqRef,ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg,needMinAmpl)
            bestMatchPrimers.append([bestMatch,theBest])
            plus40.append(seqRef[pos+1:min(len(seqRef),pos+primerLen+primerLenDev+1)])
        else: # if strand is minus
            bestPrimers.append([revComplement(seqRef[pos+2:min(len(seqRef),pos+primerLen+primerLenDev)])+pe[0],-1])
            wtPrimers.append(revComplement(seqRef[pos+2:min(len(seqRef),pos+primerLen+primerLenDev)])+wtpe[0])
            bestMatch,theBest=constructPrimers(bestPrimers[-1],wtPrimers[-1],revComplement(seqRef),ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg,needMinAmpl)
            bestMatchPrimers.append([bestMatch,theBest])
            plus40.append(revComplement(seqRef[max(pos-primerLen-primerLenDev,0):pos]))
    results=[]
    for i,bp in enumerate(bestPrimers):
        if bp[1]>0:
            if bestMatchPrimers[i][0]:
                results.append([seqNames[-1],ref,alt,'+']+bestMatchPrimers[i][0][0].split('_')+[str(primerVals[i])]+bestMatchPrimers[i][0][2])
##                rFile.write(seqNames[-1]+'\t'+ref+'/'+alt+'\t'+'+\t'+'\t'.join(bestMatchPrimers[i][0][0].split('_'))+'\t'+str(primerVals[i])+'\t'+'\t'.join(bestMatchPrimers[i][0][2])+'\n')
            else:
                results.append([seqNames[-1],ref,alt,'+']+bestMatchPrimers[i][1][0].split('_')+[str(primerVals[i])]+bestMatchPrimers[i][1][2])
##                rFile.write(seqNames[-1]+'\t'+ref+'/'+alt+'\t'+'+\t'+'\t'.join(bestMatchPrimers[i][1][0].split('_'))+'\t'+str(primerVals[i])+'\t'+'\t'.join(bestMatchPrimers[i][1][2])+'\n')
        else:
            if bestMatchPrimers[i][0]:
                results.append([seqNames[-1],revComplement(ref),revComplement(alt),'-']+bestMatchPrimers[i][0][0].split('_')+[str(primerVals[i])]+bestMatchPrimers[i][0][2])
##                rFile.write(seqNames[-1]+'\t'+revComplement(ref)+'/'+revComplement(alt)+'\t'+'-\t'+'\t'.join(bestMatchPrimers[i][0][0].split('_'))+'\t'+str(primerVals[i])+'\t'+'\t'.join(bestMatchPrimers[i][0][2])+'\n')
            else:
                results.append([seqNames[-1],revComplement(ref),revComplement(alt),'-']+bestMatchPrimers[i][1][0].split('_')+[str(primerVals[i])]+bestMatchPrimers[i][1][2])
##                rFile.write(seqNames[-1]+'\t'+revComplement(ref)+'/'+revComplement(alt)+'\t'+'-\t'+'\t'.join(bestMatchPrimers[i][1][0].split('_'))+'\t'+str(primerVals[i])+'\t'+'\t'.join(bestMatchPrimers[i][1][2])+'\n')
    return(results)
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

# Section of input arguments
par=argparse.ArgumentParser(description='endChooser automatically constructs MAMA-primers')
par.add_argument('--fasta-file','-fa',dest='fastaFile',type=str,help='Fasta-file that contains sequence with variable position designated like [A/G] or (A/G)',required=True)
par.add_argument('--amplicon-length','-alen',dest='ampliconLen',type=int,help='Amplicon length (Default: 150 bp)',required=False,default=150)
par.add_argument('--amplicon-length-deviation','-alendev',dest='ampliconLenDev',type=int,help='Amplicon length deviation (Default: 10 bp)',required=False,default=10)
par.add_argument('--primer-length','-plen',dest='primerLen',type=int,help='Optimal primer length (Default: 20 bp)',required=False,default=20)
par.add_argument('--primer-length-deviation','-plendev',dest='primerLenDev',type=int,help='Optimal primer length deviation (Default: 4 bp)',required=False,default=4)
par.add_argument('--melting-temperature','-mtemp',dest='meltTemp',type=int,help='Optimal melting temperature (Default: 60 degrees Celsius)',required=False,default=60)
par.add_argument('--melting-temperature-deviation','-mtempdev',dest='meltTempDev',type=int,help='Optimal melting temperature deviation (Default: 5 degrees Celsius)',required=False,default=5)
par.add_argument('--min-dg','-dg',dest='dimerdg',type=int,help='Minimal dG of dimer and hairpin formation (Default: 3000 kcal/mol)',required=False,default=3000)
par.add_argument('--min-amplicon','-minAmpl',dest='needMinAmpl',action='store_true',help='use this parameter if you need amplicons with a minimal length')
par.add_argument('--reference-file','-ref',dest='refFile',type=str,help='Fasta-file with the human reference genome sequence ucsc.hg19.fa. Use this parameter only if you want to check primers for covering SNPs from dbSNP',required=False)
args=par.parse_args()

fastaFile=args.fastaFile
ampliconLen=args.ampliconLen
ampliconLenDev=args.ampliconLenDev
primerLen=args.primerLen
primerLenDev=args.primerLenDev
meltTemp=args.meltTemp
meltTempDev=args.meltTempDev
dimerdg=args.dimerdg
needMinAmpl=args.needMinAmpl
refFile=args.refFile

try:
    f=open(fastaFile)
except FileNotFoundError:
    print('#######\nERROR! File not found:',fastaFile,'\n#######')
    exit(0)

print('Reading input file...')
rFile=open(fastaFile+'.primers.xls','w')
faFile=open(fastaFile+'.for_blast.fa','w')
rFile.write('Sequence_Name\tRef/Alt\tStrand\tBest_Primer_for_Mutant\tBest_Primer_for_WT\tBest_Reverse_Primer'
            '\tDiscrimination_Value\tMutant_Primer_Len\tWT_Primer_Len\tReverse_Primer_Len\tAmplicon_Length'
            '\tMutant_Primer_Tm\tWT_Primer_Tm\tReverse_Primer_Tm'
            '\tMutant_Primer_Homodimer\tWT_Primer_Homodimer\tReverse_Primer_Homodimer'
            '\tMutant_Primer_Hairpin1\tWT_Primer_Hairpin2\tReverse_Primer_Hairpin'
            '\tMutant_Rev_Primers_Heterodimer\tWT_Rev_Primers_Heterodimer\n')
seqNames=[]
print('Constructing primers...')
p=re.compile('([\(\[](\w)\/\w[\)\]])')
results=[]
for line in f:
    if '>' in line:
        if len(seqNames)>0:
            results.extend(chooseBestPrimers(seq,seqNames,rFile,ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg,needMinAmpl))
            m=p.findall(seq)
            newSeq=seq.replace(m[0][0],m[0][1])
            faFile.write(newSeq+'\n')
        seqNames.append(line.replace('> ','').replace('>','').replace('\n','').replace('\r',''))
        faFile.write(line)
        seq=''
    else:
        seq+=line.replace('\n','').replace('\r','').replace(' ','').replace('\t','').upper()
results.extend(chooseBestPrimers(seq,seqNames,rFile,ampliconLen,ampliconLenDev,primerLen,primerLenDev,meltTemp,meltTempDev,dimerdg,needMinAmpl))
m=p.findall(seq)
newSeq=seq.replace(m[0][0],m[0][1])
faFile.write(newSeq+'\n')
faFile.close()
if refFile:
    print('Checking primer sequences for covering SNPs from dbSNP...')
    blastResultFileName=faFile.name[:-2]+'blast_result.xls'
    seqsCoords={}
    cmdResult=sp.check_output('blastn -query '+faFile.name+' -subject '+refFile+' -outfmt "6 qseqid sallseqid sstrand qstart qend sstart send pident qseq sseq qlen" -perc_identity 95 -qcov_hsp_perc 95 > '+blastResultFileName,shell=True)
    file=open(blastResultFileName)
    for string in file:
        cols=string.replace('\n','').split('\t')
        chrom=cols[1]
        start=cols[5]
        end=cols[6]
        if cols[0] not in seqsFound.keys():
            seqsCoords[cols[0]]=[[chrom,start,end,qstart,qend]]
        else:
            seqsCoords[cols[0]].append([chrom,start,end,qstart,qend])
    if len(seqsCoords.keys())==0:
        print('ERROR! No sequences were found in the reference sequence:',refFile)
        exit(0)
    for i,seqName in enumerate(seqNames):
        if seqName not in seqsCoords.keys():
            print('WARNING! No sequences were found in the reference sequence',refFile,'for the sequence',seqName)
        elif len(seqsCoords[seqName])>1:
            print('WARNING! Input sequence',seqName,'has repeats in the reference sequence!')
        else:
            for res in results[2*i:2*i+2]:
                for primerStart in res[-3:]:
                    print(primerStart)
                input()
rFile.close()
print('Done')
