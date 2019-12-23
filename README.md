# spacer blasting
This tutorial will show you how to extract spacers from reads which needs UNIX, macOS, or bash on Windows.
## Required software
BLAST+

python2

## Preparing your dataset
After sequencing some genomes you'll have some paired end reads. You will trim the fastq files any way you want. You can use prinseq, trimmomatic, etc. 

trimmomatic -threads 4 -phred33 -trimlog trim.log sample_R1.fastq sample_R2.fastq sample_R1_trimmed.fastq sample_R1_trimmed_unpaired.fastq sample_R2_trimmed.fastq sample_R2_trimmed_unpaired.fastq ILLUMINACLIP:$adapters:$trimparams 2>sample_run_trimmomatic.log && rm trim.log;

MISC_TrimParameters     2:36:10 LEADING:20 TRAILING:20 MAXINFO:50:0.97 MINLEN:50


After this, you will have to convert them into fasta files from fastq files. 
You can using anything to do this such as the fastx toolkit. Or you can just use this command:

`sed -n '1~4s/^@/>/p;2~4p' fastqfile > fastafile`

if you have a paired end reads separated as to fastas, just concatenate those files into a single fasta file. And if you have spacers in your headers they should be removed to make a blastdb. To do this just type:

`sed 's/ /_/g' -i fastafile`

Now you gotta make a blastdatabase. This may take a few minutes.

`makeblastdb -in fastafile -dbtype nucl -parse_seqids`

## Extracting spacers
All you need are fasta files of repeats of interest and this blast database and we can start extracting.
There are 4 arguments you need. -t is the prefix of your spacer headers. -d is the blastdatabase you created. -q is the fasta for repeats. -o is the output.

`spacer2blast.py -t taxon -q repeat.fasta -d blastdb -o outputfile`
