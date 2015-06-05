#!/usr/bin/python
import os
import sys
import glob

def readframe(aa, pos):
  wtaa  = aa[0]
  mutaa = aa[1]
  pos   = (pos-11)/3+1
  return wtaa+str(pos)+mutaa

def adjustmutpos(muts, offset):
  if muts == 'WT': return muts
  else:
    muts  = muts.rsplit('-')
    mlist = []
    for m in muts:
      pos = str(int(m[1:-1])+offset-1)
      mlist.append(m[0]+pos+m[-1])
    return '-'.join(mlist)

def list2string(l):
  newl = []
  for i in l:
    newl.append(str(i))
  return newl

def adjustbgmuts(muts,amp,bgmuts):
  if amp in bgmuts.keys():
    for m in bgmuts[amp]:
      wtaa= m[0]
      pos  = m[1:-1]
      mutaa = m[-1]
      if muts == 'WT':
        muts = muts.replace('WT',mutaa+pos+wtaa)
      elif wtaa+pos not in muts:
        muts = muts+'-'+mutaa+pos+wtaa
      elif m in muts:
        if '-' in muts:
          muts = muts.replace('-'+m,'')
          muts = muts.replace(m+'-','')
        else:
          muts = muts.replace(m,'WT')
      elif wtaa+pos in muts:
        muts = muts.replace(wtaa+pos,mutaa+pos)
      else:
        print 'Other Conditions Exist!!!!???'
      return muts
  else: 
    return muts
      

#READ IN OFFSET FILE
infile  = open('Fasta/flu3offset','r')
offsetH = {}
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  offsetH[line[0]] = [int(line[1]),int(line[2])]
infile.close()
#READ IN BARCODE FILE
infile   = open('Fasta/BarCode','r')
barcodes = {}
pops     = []
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  barcodes[line[0]] = line[1]
  pops.append(line[1])
infile.close()

#GENOTYPE AND DEPTH HASH INITIATE
GenotypeH = {}
depthH    = {}
WTcount   = {}
for pop in pops:
  GenotypeH[pop] = {}
  depthH[pop]    = {}
  WTcount[pop]   = {}
  for amp in offsetH.keys():
    depthH[pop][amp] = 0
    WTcount[pop][amp] = 0

##############################MAIN##################################
infile = open('result/AllM_flu3','r')
bgmuts = {'flu3amp3':['C174G']}#INDICATE ALL BACKGROUND MUTATIONS
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  if line[1] not in barcodes.keys(): continue
  amp    = line[0]
  offset = offsetH[amp][0]
  pop    = barcodes[line[1]]
  muts   = line[2]
  #muts   = adjustbgmuts(muts,amp,bgmuts)
  muts   = adjustmutpos(muts,offset)
  if 'C682G' in muts: print muts, t
  depthH[pop][amp] += 1
  if muts == 'WT':
    WTcount[pop][amp] += 1
  else:
    if GenotypeH[pop].has_key(muts): 
      GenotypeH[pop][muts] += 1
    else: 
      GenotypeH[pop][muts] = 1
infile.close()

#SUMMARIZE GENOTYPE
Genotypes = []
for pop in pops:
  Genotypes.extend(GenotypeH[pop].keys())
Genotypes = list(set(Genotypes))

#OUTPUT
#WTCOUNT AND DEPTH OUTPUT
outfile = open('result/WT_N_Depth','w')
outfile.write('pop'+"\t"+'amp'+"\t"+'WTcount'+"\t"+'Depth'+"\n")
for pop in pops:
  for amp in offsetH.keys():
    outfile.write("\t".join([pop, amp, str(WTcount[pop][amp]), str(depthH[pop][amp])])+"\n")
outfile.close()

#GENOTYPES
outfileA = open('result/AGenotypes','w')
outfileS = open('result/SGenotypes','w')
header   = 'Genotype'+"\t"+"\t".join(pops)
outfileA.write(header+"\n")
outfileS.write(header+"\n")
for G in Genotypes:
  counts = []
  for pop in pops:
    if GenotypeH[pop].has_key(G): 
      counts.append(GenotypeH[pop][G])
    else: 
      counts.append(0)
  out = G+"\t"+"\t".join(list2string(counts))
  outfileA.write(out+"\n")
  if '-' not in G:
    outfileS.write(out+"\n")
outfileA.close()
outfileS.close()

#SINGLE MUTS INFO FILE
outfileC = open('result/SMut','w')
header   = ['Pos','Genotype','Frame1','Frame2']
for pop in pops:
  header.extend([pop, pop+'_WT',pop+'_Dep'])
outfileC.write("\t".join(header)+"\n")

#READ IN PROTEIN LEVEL DATA
infile = open('Fasta/flu3info','r')
Rframe = {}
for line in infile.xreadlines():
  array = line.rstrip().rsplit("\t")
  ID = ''.join(array[0:3])
  WTaa1 = array[3]
  Mutaa1 = array[4]
  WTaa2 = array[5]
  Mutaa2 = array[6]
  Rframe[ID] = [WTaa1+Mutaa1, WTaa2+Mutaa2]
infile.close()

for G in Genotypes:
  if '-' not in G and 'N' not in G:
    pos = int(G[1:-1])
    F1  = readframe(Rframe[G][0],pos)
    F2  = readframe(Rframe[G][1],8)
    out = [str(pos),G,F1,F2]
    for pop in pops:
      wtc = 0
      dep = 0 
      for amp in offsetH.keys():
        if pos >= offsetH[amp][0] and pos <= offsetH[amp][1]:
          wtc += WTcount[pop][amp]
          dep += depthH[pop][amp]
      if GenotypeH[pop].has_key(G): out.append(str(GenotypeH[pop][G]))
      else: out.append(str(0))
      out.append(str(wtc))
      out.append(str(dep))
    outfileC.write("\t".join(out)+"\n")
outfileC.close()
