# spacer blasting
This tutorial will show you how to extract spacers from reads.
##Required software
BLAST+

## Preparing your dataset
After sequencing some genomes you'll have some paired end reads. You will trim the fastq files any way you want. You can use prinseq, trimmomatic, etc. After this, you will have to convert them into fasta files from fastq files. 
You can using anything to do this such as the fastx toolkit. Or you can just use this command:

`sed -n '1~4s/^@/>/p;2~4p' fastqfile > fastafile`

if you have a paired end reads separated as to fastas, just concatenate those files into a single fasta file. Now you gotta make a blastdatabase. This may take a few minutes.

`makeblastdb -in fastafile -dbtype nucl -parse_seqids`

## Extracting spacers
All you need are fasta files of repeats of interest and this blast database and we can start extracting.
There are 4 arguments you need. -t is the prefix of your spacer headers. -d is the blastdatabase you created. -q is the fasta for repeats. -o is the output.

`spacer2blast.py -t taxon -q repeat.fasta -d blastdb -o outputfile`
