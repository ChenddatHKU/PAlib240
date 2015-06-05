#!/usr/bin/python
import os
import sys
import glob
from string import atof

def floor(n):
  if atof(n) < 0.001: return str(0.001)
  else: return n

infile  = open('result/SMut','r')
outfile = open('result/HCMut','w')
silfile = open('tmp/sil','w')
misfile = open('tmp/mis','w')
nonfile = open('tmp/non','w')
header  = "\t".join(['Mut','AA','DNA_A','DNA_B','Tra_A','Tra_B','Inf_A1','Inf_A2','Fit'])
outfile.write(header+"\n")
countline  = 0
mincov     = 15000
count      = 0
counthc    = 0
for line in infile.xreadlines():
  if countline == 0: countline += 1; continue
  line       = line.rstrip().rsplit("\t")
  Pos        = int(line[0])
  G          = line[1]
  F1         = line[2]
  F2         = line[3]
  WT         = atof(line[4])
  WT_WT      = atof(line[5])
  WT_Dep     = atof(line[6])
  DNA_A      = atof(line[7])
  DNA_A_WT   = atof(line[8])
  DNA_A_Dep  = atof(line[9])
  DNA_B      = atof(line[10])
  DNA_B_WT   = atof(line[11])
  DNA_B_Dep  = atof(line[12])
  Tra_A      = atof(line[13])
  Tra_A_WT   = atof(line[14])
  Tra_A_Dep  = atof(line[15])
  Tra_B      = atof(line[16])
  Tra_B_WT   = atof(line[17])
  Tra_B_Dep  = atof(line[18])
  Inf_A1     = atof(line[19])
  Inf_A1_WT  = atof(line[20])
  Inf_A1_Dep = atof(line[21])
  Inf_A2     = atof(line[22])
  Inf_A2_WT  = atof(line[23])
  Inf_A2_Dep = atof(line[24])
  if DNA_A_Dep < mincov or DNA_B_Dep < mincov or Tra_A_Dep < mincov or Tra_B_Dep < mincov or Inf_A1_Dep < mincov or Inf_A2_Dep < mincov: continue
  WTfreq  = WT/WT_Dep
  DNAfreq = (DNA_A+DNA_B)/(DNA_A_Dep+DNA_B_Dep)
  Trafreq = Tra_A/Tra_A_Dep
  count += 1
  if DNAfreq < WTfreq*4: continue
  if DNAfreq < 5/atof(mincov): continue
  counthc += 1
  DNA_A_RF  = str(DNA_A/DNA_A_WT)
  DNA_B_RF  = str(DNA_B/DNA_B_WT)
  Tra_A_RF  = str(Tra_A/Tra_A_WT)
  Tra_B_RF  = str(Tra_B/Tra_B_WT)
  Inf_A1_RF = str(Inf_A1/Inf_A1_WT)
  Inf_A2_RF = str(Inf_A2/Inf_A2_WT)
  DNA_RF    = (DNA_A+DNA_B)/(DNA_A_WT+DNA_B_WT)
  Tra_RF    = Tra_A/Tra_A_WT
  Inf_RF    = (Inf_A1+Inf_A2)/(Inf_A1_WT+Inf_A2_WT)
  fit = str((atof(Inf_A1_RF)+atof(Inf_A2_RF))/DNA_RF/2)
  out = "\t".join([G,F1,DNA_A_RF,DNA_B_RF,Tra_A_RF,Tra_B_RF,Inf_A1_RF,Inf_A2_RF,fit])
  outfile.write(out+"\n")
  if Pos > 400 and Pos < 1800:
    if F1[-1] == '_': nonfile.write(str(floor(fit))+"\n")
    elif F1[0] != F1[-1]: misfile.write(str(floor(fit))+"\n")
    elif F1[0] == F1[-1]: silfile.write(str(floor(fit))+"\n")
infile.close()
silfile.close()
misfile.close()
nonfile.close()
#print count, counthc

