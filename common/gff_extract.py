#Extract data required for gene plots from a gff file (from a single gene gff)
import os

x=[]
exon_start=''
exon_end=''
numexon=0

a=open("/Users/shrirambhat/Desktop/telomere/data/briggsae/gffs/genes/qx_9296.gff3","r")
for i in a: 
	if i[0]=='#':
		continue
	x.append(i.split('\t')) #x is temporary

chr=x[0][0] 

for i in x:
	for j in range(len(i)):
		if i[j]=='exon':
			exon_start = exon_start + i[j+1] + ','
			exon_end = exon_end + i[j+2] + ','
			numexon+=1

exon_start=exon_start.strip(',')
exon_end=exon_end.strip(',')


print('exon start: '+exon_start)
print('exon end: '+exon_end)
print('number of exons: '+str(numexon))
print('chromsome number: '+chr)