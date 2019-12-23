#!/usr/bin/python

import sys
import getopt
import optparse
from optparse import OptionParser
from collections import defaultdict
import subprocess
import os

def hamming2(str1, str2):
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
           diffs += 1
    return diffs


def makeRC(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a':'t','c':'g','g':'c','t':'a'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def callblastn(db, query):
    outfile=query+'_vs_'+os.path.basename(db)
    subprocess.call(['blastn', '-db', db , '-task', 'blastn-short', '-query', query, '-num_threads','6', '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop qlen slen', '-out', outfile, '-max_target_seqs', '1000000', '-evalue', '0.001'])
    return outfile 
def parsespacers(blastfile):
    with open(blastfile, "r") as br:
        query=''
        subject2startend={}
        qsubdic={}
        for line in br:
            linearray=line.split('\t')
            if linearray[0] != query:
                subject2startend={}
                query=linearray[0]
            subject=linearray[1]
            sstart,send=int(linearray[8]),int(linearray[9])
            subject2startend.setdefault(subject,[]).append([sstart,send])
            qsubdic[query]=subject2startend
        return qsubdic
                
def sortdiccoords(qsubdic,db):
    batchlist=[]
    for repeats in qsubdic:
        fname=repeats+'.'+os.path.basename(db)+'.entry'
        batchlist.append(fname)
        f=open(fname, "w")
        for reads in qsubdic[repeats]:
            qsubdic[repeats][reads].sort(key=lambda x:x[0])
            if len(qsubdic[repeats][reads])>1:
                #print repeats, reads, qsubdic[repeats][reads]
                locs=qsubdic[repeats][reads]
                x=0
                while len(locs)-1>x:
                    s1,s2,e1,e2=locs[x][0],locs[x+1][0],locs[x][1],locs[x+1][1]
                    if s1<e1 and s2-e1<50 and s2-e1>20:
                        line=' '.join([reads, str(e1+1)+'-'+str(s2-1),'plus','\n'])
                        f.write(line)
                    if s1>e1 and e2-s1<50 and e2-s1>20:
                        line=' '.join([reads, str(s1+1)+'-'+str(e2-1),'minus','\n'])
                        f.write(line)
                    x=x+1
        f.close()
        if os.stat(fname).st_size==0:
            os.remove(fname)
            batchlist.remove(fname)
    return batchlist

def createspacerarray(fastafile,taxon):         
    taxon2spacer={}
    with open(fastafile, 'r') as spacerfile:
        for line1 in spacerfile:
            if '>' not in line1:
                taxon2spacer.setdefault(taxon, []).append(line1.strip())
    return taxon2spacer
                
def blastdbcmd(entrybatch, db):
    spacerfasta=entrybatch+'.fasta'
    subprocess.call(['blastdbcmd', '-entry_batch', entrybatch, '-db', db, '-out', spacerfasta])
    return spacerfasta
    
#def tempfasta(spacerlist):
#    f=open('tempfasta', 'w')
#    for i in range(len(set(spacerlist))):
#        f.write('>'+str(i)+'\n')

def renamespacerfile(taxon2spacer,filename):
    spacerstatfname='spacerstats.csv'
    stat=open(spacerstatfname,"w")
    for taxon in taxon2spacer:
	spacerlist=list(set(taxon2spacer[taxon]))
	print taxon, len(spacerlist)
	counterlist=taxon2spacer[taxon]
	for spacers in spacerlist:
	    for spacers2 in spacerlist:
		if spacers != spacers2:
		    if hamming2(spacers, spacers2)/float(len(spacers)) < .1 or hamming2(spacers,makeRC(spacers2))/float(len(spacers))<.1:
			if counterlist.count(spacers)>counterlist.count(spacers2):
			    spacerlist=filter(lambda a: a !=spacers2, spacerlist)
			   # spacerlist.remove(spacers2)
			elif counterlist.count(spacers)<counterlist.count(spacers2):
			    spacerlist=filter(lambda a: a !=spacers, spacerlist)
			    #spacerlist.remove(spacers)
			else:
			    spacerlist=filter(lambda a: a !=spacers2, spacerlist)
		    if hamming2(spacers, spacers2[1:])/float(len(spacers)) < .1 or hamming2(spacers,makeRC(spacers2[1:]))/float(len(spacers))<.1:
			if counterlist.count(spacers)>counterlist.count(spacers2):
			    spacerlist=filter(lambda a: a !=spacers2, spacerlist)
			   # spacerlist.remove(spacers2)
			elif counterlist.count(spacers)<counterlist.count(spacers2):
			    spacerlist=filter(lambda a: a !=spacers, spacerlist)
			    #spacerlist.remove(spacers)
			else:
			    spacerlist=filter(lambda a: a !=spacers2, spacerlist)
#		elif spacers == spacers2:
#			    spacerlist=filter(lambda a: a !=spacers2, spacerlist)
			    #spacerlist.remove(spacers2)
	#taxon2spacer[taxon]=spacerlist
	statstring=','.join([taxon, str(len(spacerlist)), str(len(taxon2spacer[taxon])),'\n']) 
	stat.write(statstring)
	print len(spacerlist)
	f=open(filename, "w")
	spacercount=1
        taxon2spacer[taxon]=spacerlist
	for finalspacer in taxon2spacer[taxon]:
	    spacerheader='>'+taxon+'_'+str(spacercount)+'\n'
	    spacercount=spacercount+1
	    f.write(spacerheader)
	    finalspacerstring=finalspacer+'\n'
	    f.write(finalspacerstring)
	#print cmd_out
    stat.close()
    f.close()
def main():
    parser=OptionParser(usage="spacer2blast.py -t taxon -q repeat.fasta -d blastdb -o outfile") 
    parser.add_option("-q", "--query", dest="q",  help="queryfile") 
    parser.add_option("-d", "--dbfile", action="store", dest="db", default=False,help="dbfile") 
    parser.add_option("-o", "--outputfile", action="store", dest='out', default=False,help="outputfile") 
    parser.add_option("-t", "--taxon", action="store", dest='t', default=False,help="taxon") 
    (options,args)=parser.parse_args() 
    if len(args) !=0: 
        parser.error("need input and output") 
    blastreport=callblastn(options.db, options.q)
    qsubdic=parsespacers(blastreport)    
    batchlist=sortdiccoords(qsubdic,options.db)
    for entrybatch in batchlist:
        spacerfasta=blastdbcmd(entrybatch,options.db)
        print spacerfasta
        tax2spacer=createspacerarray(spacerfasta,options.t)
        renamespacerfile(tax2spacer, options.out)
    if len(args) !=0: 
        parser.error("need input and output") 
    #callblastn(db, queryfile,out)
    

#callblastn(sys.argv[1], sys.argv[2], sys.argv[3]

if __name__== '__main__':
    main()
    #callblastn(db,options.out)

