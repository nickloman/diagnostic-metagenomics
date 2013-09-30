read1=$1
read2=$2
out=$3

bin/bowtie2-2.0.0-beta5/bowtie2 --very-fast -p 8 -x refs/human/hg19 -1 $read1 -2 $read2 | \
	bin/samtools-0.1.18/samtools view -f 4 -S -h - | \
	bin/samtools-0.1.18/samtools view -f 8 -S -h -b -o $3 -


