import sys
import subprocess

def callblastn(subject):
    #outfile=query+'_vs_'+os.path.basename(db)
    outfile='temp.aln'
    subprocess.call(['blastn', '-subject', subject , '-task', 'blastn-short', '-query', subject, '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop qlen slen', '-out', outfile, '-max_target_seqs', '1000000', '-evalue', '0.001'])
    return outfile 

callblastn(sys.argv[1])

header2fasta={}
headerarray=[]
with open(sys.argv[1]) as fasta:
    for line in fasta:
        if '>' in line:
            header=line.strip()[1:]
            headerarray.append(line.strip()[1:])
    	else:
            header2fasta[header]=line.strip()

headerarray2=[]
with open('temp.aln') as alnfile:
    for line in alnfile:
      if int(line.split('\t')[3])> 20:
 #       if line.split('\t')[0]==line.split('\t')[1]:w
        if line.split('\t')[0] != line.split('\t')[1]:
          if line.split('\t')[0] in headerarray:
            if line.split('\t')[0] not in headerarray2:
               headerarray.remove(line.split('\t')[0])
               headerarray2.append(line.split('\t')[1])

for header in headerarray:
    print '>'+header
    print header2fasta[header]
        
	    
