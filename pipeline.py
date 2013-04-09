# export PATH=$PATH:/home/nick/ecoli_metagenomics/g:/mnt/fast/software/aligners/samtools-0.1.17/:/home/nick/ecoli_metagenomics/bin/bowtie2-2.0.0-beta5/


import os
import re
from ruffus import *
from runutils import read_run_details
import subprocess
import sqlite3

parser = cmdline.get_argparse(description='Clinical metagenomics pipeline')

parser.add_argument("-s", "--sqlite_db", type=str, help="Name and path of SQLite3 database")
parser.add_argument("-c", "--sql_command", type=str, help="SQL command to return rows")

options = parser.parse_args()

#  optional logger which can be passed to ruffus tasks
logger, logger_mutex = cmdline.setup_logging (__name__, options.log_file, options.verbose)

#_____________________________________________________________________________________
#   pipelined functions go here
#_____________________________________________________________________________________

def cmd_exists(cmd):
    try:
        ret = subprocess.call([cmd], 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return ret == 0
    except OSError, e:
        return 0

if not cmd_exists("bowtie2"):
	raise SystemExit("Bowtie2 not on path")

# runs = read_run_details("runs.txt")

def read_run_details_sqlite3(fn, cmd):
	conn = sqlite3.connect(fn)
	conn.row_factory=sqlite3.Row
	c = conn.cursor()
	c.execute(cmd)
	rows = c.fetchall()
	return rows

runs = read_run_details_sqlite3(options.sqlite_db, options.sql_command)
print "%d runs selected " % (len(runs),)

usefulhash = {}

initial_input_files = []
for run in runs:
    infiles = []
    for n, read in enumerate([str(x) for x in range(1,29)]):
        readid = 'File' + read
        if len(run[readid]):
            path = '%s/%s' % (run['RunFolder'], run[readid])
            infiles.append(path)

    out = '%s.fastq.gz' % (run['Label'])
    label = run['Label']
    short_label = run['ShortLabel']
    initial_input_files.append([infiles, out, label, short_label])
    rundict = dict(run)
    rundict['label'] = label
    usefulhash[out] = rundict
    usefulhash[label] = rundict

def get_labels(infiles):
    labels = []
    for f in infiles:
       for f2 in initial_input_files:
           if f.startswith(f2[2]):
               labels.append((f2[3], f))
    if len(labels) != len(infiles):
       raise SystemExit("SILLY!")
    return labels

# @files('qc.R', 'cov.pdf')

@files(initial_input_files)
def merge_fastq(infiles, outfiles, labels, short_labels):
    cmd = "zcat %s | gzip > %s" % (" ".join(infiles), outfiles)
    os.system(cmd)

@merge(merge_fastq, 'results/qc.pdf')
def qc_on_reads(infiles, outfile):
    labels = [usefulhash[f]['label'] for f in infiles]
    files = zip(infiles, labels)
    script = ("Rscript qc.R %s" % (" ".join( ["%s=%s" % (i[1], i[0]) for i in files])))
    print script
    os.system(script)

@jobs_limit(1)
@transform(merge_fastq, suffix(".fastq.gz"), ".bam")
def filter_out_human_reads(infile, outfile):
    os.system("scripts/filter_human_reads_fragment.sh %s %s %s.stats.txt" % (infile, outfile, usefulhash[infile]['label']))

@merge(filter_out_human_reads, "results/outstats.txt")
def collect_human_match_stats(infiles, outfile):
    os.system("python scripts/collect_match_stats.py %s" % (options.sqlite_db))

@transform(filter_out_human_reads, suffix(".bam"), ".fastq")
def bam_to_fastq(infile, outfile):
    os.system("bin/bam2fastq-1.1.0/bam2fastq -o %s %s" % (outfile, infile))
   
@jobs_limit(1) 
@transform(bam_to_fastq, suffix(".fastq"), ".fastq.bowtie2out.txt")
def taxonomic_assignment_by_clade_specific_genes(infile, outfile):
    os.system("bin/nsegata-metaphlan-8485393d6b53/metaphlan.py %s --bowtie2db bin/nsegata-metaphlan-8485393d6b53/bowtie2db/mpa --nproc 8 > results/%s.relreport.txt" % (infile, infile))

@jobs_limit(1) 
@transform(bam_to_fastq, suffix(".fastq"), ".fastq.bowtie2out.176.txt")
def taxonomic_assignment_by_clade_specific_genes_metaphlan_176(infile, outfile):
    os.system("bin/metaphlan-1.7.6/metaphlan.py %s --bowtie2db bin/metaphlan-1.7.6/mpa --nproc 8 --bowtie2out %s > results/%s.relreport-176.txt" % (infile, outfile, infile))


@transform(taxonomic_assignment_by_clade_specific_genes, suffix(".fastq.bowtie2out.txt"), ".cladeprofile.txt")
def metaphlan_cladeprofile(infile, outfile):
    os.system("bin/nsegata-metaphlan-8485393d6b53/metaphlan.py -t clade_profiles %s > results/%s.cladeprofile.txt" % (infile, infile))

@jobs_limit(1)
@transform(bam_to_fastq, suffix(".fastq"), ".soapdenovo.contig.stats.txt")
def soapdenovo(infile, outfile):
    tag = infile[0:-6]
    cfg = """
max_rd_len=151
[LIB]
reverse_seq=0
asm_flag=1
q=%s""" % (infile,)
    cfg_file = ("%s.soapdenovo.txt" % (tag,))
    with open(cfg_file, "w") as f:
        f.write(cfg)
    os.system("./soapgo %s %s.soapdenovo 39" % (infile, tag))

@jobs_limit(1)
@merge(bam_to_fastq, "coassembly/coassembly.fastq")
def collate_for_assembly(infiles, outfile):
    os.system("cat %s > %s" % (" ".join(infiles), outfile))

#@jobs_limit(1)
#@transform(bam_to_fastq, suffix(".fastq"), "_vs_280.bam")
#def map_to_stec(infile, outfile):
#    os.system("bwa bwasw -t 8 refs/280_flxplus.fna %s | samtools view -F 4 -o %s -S -b -" % (infile, outfile))

@jobs_limit(1)
@transform(bam_to_fastq, suffix(".fastq"), r'aligned_2011c/\1_vs_2011c.sorted.bam')
def map_to_stec(infile, outfile):
    os.system("bwa bwasw -t 8 refs/2011c.fna %s | samtools view -F 4 -S -b - | samtools sort - %s" % (infile, outfile[0:-4]))

@jobs_limit(1)
@transform(bam_to_fastq, suffix(".fastq"), r'coassembly/aligned/\1_vs_coassembly.sorted.bam')
def map_to_coassembly(infile, outfile):
    os.system("bwa bwasw -t 8 coassembly/coassembly_10pc.ray/Contigs.fasta %s | samtools view -F 4 -S -b - | samtools sort - %s" % (infile, outfile[0:-4]))

@jobs_limit(1)
@transform(bam_to_fastq, suffix(".fastq"), r'coassembly/aligned_steccontigs/\1_vs_steccontigs.sorted')
def map_to_steccontigs(infile, outfile):
    os.system("bwa bwasw -t 16 coassembly/cluster10pc/stec20contigs.fasta %s | samtools view -F 4 -S -b - | bamtools sort -in /dev/stdin -out %s" % (infile, outfile))

@transform(map_to_stec, suffix(".bam"), ".sorted.bam")
def sort_stec(infile, outfile):
    os.system("samtools sort %s %s.sorted" % (infile, infile[0:-4]))

@transform(sort_stec, suffix(".sorted.bam"), ".depths.Rdata")
def depths(infile, outfile):
    os.system("Rscript scripts/dumpcov.R %s scaffold00001 %s" % (infile, outfile))
    #os.system("samtools depth %s -Q30 > %s.Q30" % (infile, outfile))
    #os.system("samtools depth %s > %s" % (infile, outfile))

@transform(sort_stec, suffix(".sorted.bam"), ".sorted.bam.bai")
def index_stec(infile, outfile):
    os.system("samtools index %s" % (infile))

@jobs_limit(1)
@transform(bam_to_fastq, suffix(".fastq"), r'aligned_hmpdacc/\1_vs_hmpdacc.bam')
def map_to_hmpdacc(infile, outfile):
    os.system("bwa bwasw -t 16 refs/hmpdacc/bact_seqs.fa %s | samtools view -o %s -S -b -" % (infile, outfile))

@transform(map_to_hmpdacc, suffix(".bam"), ".last.txt")
def last_against_nr(infile, outfile):
    os.system("samtools view -f 4 %s | cut -f 1,10 | xargs -L 2 printf \">%%s\\n%%s\\n\" | parallel --pipe --recstart '>' /mnt/twinterry/ecoli_metagenomics/bin/last-278/src/lastal -F15 -f 0 -v /mnt/phatso/blastdb/nr.last - > %s" % (infile, outfile))

@transform(bam_to_fastq, suffix(".fastq"), ".fasta.gz")
def fastq_to_trimmed_fasta(infile, outfile):
    os.system("bin/lh3-seqtk-771d60b/seqtk trimfq %s | bin/lh3-seqtk-771d60b/seqtk seq -A /dev/stdin | gzip > %s" % (infile, outfile))

@transform(fastq_to_trimmed_fasta, suffix(".fasta.gz"), ".silva.blast.txt.gz")
def taxonomic_assignment_to_silva(infile, outfile):
    os.system("bwa bwasw refs/silva/SSURef_108_tax_silva.fasta %s | samtools view -S -F 4 - | cut -f 1 | python scripts/spitout.py %s | /mnt/fast/software/blast/bin/blastall -p blastn -F \"m D\" -m 8 -d refs/silva/SSURef_108_tax_silva.fasta -i /dev/stdin -a 1 | gzip > %s" % (infile, infile, outfile)) 

@jobs_limit(1)
@transform(fastq_to_trimmed_fasta, suffix(".fasta.gz"), ".blastx.nr.txt")
def blastx_nr(infile, outfile):
    os.system("zcat %s | /mnt/fast/software/blast-2.2.22/bin/blastall -p blastx -i /dev/stdin -d ../blastdb/nr -m 8 -a 7  > %s" % (infile, outfile))

@transform(fastq_to_trimmed_fasta, suffix(".fasta.gz"), ".silva.blat.txt.gz")
def taxonomic_assignment_to_silva_blat(infile, outfile):
    os.system("bwa bwasw refs/silva/SSURef_111_NR_tax_silva_trunc.fasta %s | samtools view -S -F 4 - | cut -f 1 | python scripts/spitout.py %s > %s.tmp" % (infile, infile, outfile))
    os.system("../bin/x86_64/blat -q=dna -t=dna -out=blast8 refs/silva/SSURef_108_tax_silva.dna.fasta %s.tmp %s" % (outfile, outfile))

#@collate(map_to_stec, regex(r'(.+)-.*'), r'\1.bam')
#def merge_bam_files(infile, outfile):
#    if len(infile) == 2:
#       os.system("samtools merge %s %s" % (outfile, " ".join(infile)))
#    else:
#       os.system("ln -s %s %s" % (infile[0], outfile))

@jobs_limit(1)
# @collate(bam_to_fastq, regex(r'(.+)-.*'), r'\1.vfdb.sam')
@transform(bam_to_fastq, suffix(".fastq"), ".vfdb.sam")
def map_reads_to_virulence_db(infile, outfile):
    os.system("cat %s | bwa bwasw -t 8 vfdb/CP_VFs.ffn /dev/stdin | samtools view -h -S -F 4 /dev/stdin > %s" % (infile, outfile))

@transform(map_to_stec, suffix(".bam"), ".covplot.pdf")
def covplot(infile, outfile):
    os.system("scripts/cov.R %s %s" % (infile, outfile))

@merge(map_reads_to_virulence_db, "results/virreport.txt")
def virulence_genes_of_interest(infile, outfile):
    os.system("python scripts/read_counts_from_sam.py virulence_genes_of_interest.txt %s > %s" % (" ".join(infile), outfile))

@jobs_limit(1)
@transform(bam_to_fastq, suffix(".fastq"), "_vs_pangenome.bam")
def map_to_pangenome(infile, outfile):
    os.system("bwa bwasw -t 8 pangenome/all_ecoli_shig_pangenome.fasta %s | samtools view -F 4 -o %s -S -b -" % (infile, outfile))

@transform(map_to_pangenome, suffix(".bam"), ".sorted.bam")
def sort_pangenome(infile, outfile):
    os.system("samtools sort %s %s.sorted" % (infile, infile[0:-4]))

@merge(sort_pangenome, 'results/pangenomedepths.txt')
def pangenome_stats(infiles, outfile):
    labels = get_labels(infiles)
    script = "hello %s" % (" ".join( ["%s=%s" % (f[0], f[1]) for f in labels]))
    print script
    #os.system(script)

@transform(bam_to_fastq, suffix(".fastq"), ".insertsize.txt")
def pairup_calculate_insert(infile, outfiles):
    tag = infile[0:-6]
    os.system("python scripts/pair_up_reads.py %s %s" % (infile, infile[:-6]))
    os.system("bwa bwasw refs/280_flxplus.fna %s_1.fastq %s_2.fastq 1>/dev/null 2>%s.insertsize.txt" % (tag, tag, tag))

#     bowtie2 -x refs/280_flxplus.fna -U 2638-N12-STEC_V_High_interleaved_pairs.fastq --very-fast-local -p 8 --reorder --mm | 
#    tee >(samtools view -S -b - > 280_flxplus.fna.2638-N12-STEC_V_High_interleaved_pairs.fastq.bowtie2.bam) |
#    tee >(samtools view -S -b - | samtools sort -m 2000000000 - 280_flxplus.fna.2638-N12-STEC_V_High_interleaved_pairs.fastq.bowtie2.sorted) |
#    sam_len_cov_gc_insert.pl -i -f refs/280_flxplus.fna -s - -out 280_flxplus.fna.2638-N12-STEC_V_High_interleaved_pairs.fastq

# bin/bowtie2-2.0.0-beta5/bowtie2-build refs/280_flxplus.fna refs/280_flxplus.fna
# bin/velvet_1.1.06/velveth velvet_2638_k12 21  -shortPaired -fastq 2638-N12-STEC_V_High_interleaved_pairs.fastq
# bin/velvet_1.1.06/velvetg velvet_2638_k12 -ins_length 335 -ins_length_sd 30 -exp_cov auto -clean yes -unused_reads yes -read_trkg no
# bin/velvet_1.1.06/velvetg velvet_2638_k12 -ins_length 335 -ins_length_sd 30 -exp_cov auto -unused_reads yes -read_trkg no

## extra commands
# bwa bwasw -t 8 refs/ecoli/Bacteria/Campylobacter_concisus_13826_uid58667/NC_009802.fna 1253-H-Cdiff.fastq  | samtools view -uS - -F 4 -q 30 | samtools sort - 1253-H-Cdiff_vs_concisus.sorted.bam
# bwa bwasw -t 8 refs/ecoli/Bacteria/Campylobacter_jejuni_NCTC_11168___ATCC_700819_uid57587/NC_002163.fna 4961-H-Campy.fastq | samtools view -uS - -F 4 | samtools sort - 4961-H-Campy_vs_cjejuni11168.fastq.sorted
# bwa bwasw -t 8 refs/ecoli/Bacteria/Clostridium_difficile_630_uid57679/NC_009089.fna 1122-H-Cdiff.fastq | samtools view -uS - -F 4 -q 30 | samtools sort - 1122-H-Cdiff_vs_cdif630.fastq.sorted
#Rscript scripts/dumpcov.R 1253-H-Cdiff_vs_concisus.sorted.bam.bam "gi|157163852|ref|NC_009802.1|" results/covdepths/1253-H-Cdiff_vs_concisus.depths.Rdata


@transform(bam_to_fastq, suffix(".fastq"), "_vs_campy.done")
def align_campy(infile, outfile):
   files = """./Campylobacter_concisus_13826_uid58667/NC_009795.fna
./Campylobacter_concisus_13826_uid58667/NC_009802.fna
./Campylobacter_concisus_13826_uid58667/NC_009796.fna
./Campylobacter_jejuni_NCTC_11168___ATCC_700819_uid57587/NC_002163.fna
./Campylobacter_hominis_ATCC_BAA_381_uid58981/NC_009714.fna
./Campylobacter_hominis_ATCC_BAA_381_uid58981/NC_009713.fna
./Campylobacter_lari_RM2100_uid58115/NC_012039.fna
./Campylobacter_lari_RM2100_uid58115/NC_012040.fna
./Campylobacter_jejuni_S3_uid159533/NC_017281.fna
./Campylobacter_jejuni_S3_uid159533/NC_017282.fna
./Campylobacter_jejuni_ICDCCJ07001_uid61249/NC_014801.fna
./Campylobacter_jejuni_ICDCCJ07001_uid61249/NC_014802.fna
./Campylobacter_jejuni_81_176_uid58503/NC_008770.fna
./Campylobacter_jejuni_81_176_uid58503/NC_008787.fna
./Campylobacter_jejuni_81_176_uid58503/NC_008790.fna
./Campylobacter_jejuni_RM1221_uid57899/NC_003912.fna
./Campylobacter_jejuni_PT14_uid176499/NC_018709.fna
./Campylobacter_jejuni_M1_uid159535/NC_017280.fna
./Campylobacter_fetus_82_40_uid58545/NC_008599.fna
./Campylobacter_curvus_525_92_uid58669/NC_009715.fna
./Campylobacter_jejuni_81116_uid58771/NC_009839.fna
./Campylobacter_jejuni_doylei_269_97_uid58671/NC_009707.fna
./Campylobacter_jejuni_IA3902_uid159531/NC_017284.fna
./Campylobacter_jejuni_IA3902_uid159531/NC_017279.fna
./Campylobacter_jejuni_NCTC_11168_BN148_uid174152/NC_018521.fna"""
   refs = ["refs/ecoli/Bacteria/" + f for f in files.split("\n")]
   tag = infile[0:-5]
   for ref in refs:
       os.system("bwa bwasw -t 8 %s %s | samtools view -uS - -F 4 | samtools sort - %s_vs_%s.fastq.sorted" % (ref, infile, tag, os.path.basename(ref)))


@transform(bam_to_fastq, suffix(".fastq"), "_vs_salmonellae_fullhit.done")
def look_for_salmonellae(infile, outfile):
   tag = outfile[0:-5]
   #os.system("bwa bwasw -t 8 refs/ecoli/Bacteria/salmonellae.fna %s | grep 150M | grep NM:i:0 > %s.sam" % (infile, tag))
   #os.system("python scripts/sam2fasta.py < %s.sam > %s.fasta" % (tag, tag))
   #os.system("/mnt/fast/software/blast-2.2.22/bin/blastall -p blastn -d /mnt/phatso/blastdb/nt -a 8 -m 7 -i %s.fasta > %s.blast.xml" % (tag, tag))
   #os.system("python scripts/blasttophit.py Salmonella < %s.blast.xml > %s.blast.fastahdrs" % (tag, tag))
   os.system("python scripts/fetchfastq.py %s.blast.fastahdrs %s > %s.blast.pairs.fasta" % (tag, infile, tag))
   os.system("/mnt/fast/software/blast-2.2.22/bin/blastall -a 8 -p blastn -i %s.blast.pairs.fasta -d /mnt/phatso/blastdb/nt -a 8 >%s.blast.pairs.blast.txt" % (tag, tag))
   os.system("touch %s.done" % (tag))

@transform(bam_to_fastq, suffix(".fastq"), "_vs_mlst.sorted.bam")
def mlst_genotype(infile, outfile):
   os.system("bwa bwasw -t 8 mlstdb/ecoli/baits.fasta %s | samtools view -bS -F 4 - | samtools sort - %s"  % (infile, outfile[0:-4]))

@transform(mlst_genotype, suffix(".sorted.bam"), ".fasta")
def mlst_genotype2(infile, outfile):
   os.system("samtools mpileup -u %s -f mlstdb/ecoli/baits.fasta |  bcftools view -cg - | vcfutils.pl vcf2fq -d 1| bin/lh3-seqtk-771d60b/seqtk seq -A - > %s.fasta" % (infile, outfile))
   os.system("blastall -p blastn -i %s.fasta -d mlstdb/ecoli/alleles.fas > %s.blast.txt" % (outfile, outfile))

@transform(bam_to_fastq, suffix(".fastq"), "_vs_antigen_toxin.sorted.bam")
def antigen_and_toxin(infile, outfile):
   os.system("bwa bwasw -t 8 refs/H-antigen-Shiga-toxin.fasta %s | samtools view -bS -F 4 - | samtools sort - %s"  % (infile, outfile[0:-4]))

#samtools mpileup -uf mlstdb/ecoli/baits.fasta 
#samtools mpileup -l cdiff.bed -u 1122-H-Cdiff_vs_cdif630.fastq.sorted.bam -f refs/ecoli/Bacteria/Clostridium_difficile_630_uid57679/NC_009089.fna | bcftools view -cg - | vcfutils.pl vcf2fq -d 1 > cdiff.txt
#samtools mpileup -u 1253-H-Cdiff_vs_concisus.sorted.bam.bam -f refs/ecoli/Bacteria/Campylobacter_concisus_13826_uid58667/NC_009802.fna  | bcftools view -cg - | vcfutils.pl vcf2fq -d 1 | bin/lh3-seqtk-771d60b/seqtk seq -A - > results/1253-concisus_consensus.fasta
#samtools mpileup -u 4961-H-Campy_vs_cjejuni11168.fastq.sorted.bam  -f refs/ecoli/Bacteria/Campylobacter_jejuni_NCTC_11168___ATCC_700819_uid57587/NC_002163.fna  | bcftools view -cg - | vcfutils.pl vcf2fq -d 1 | bin/lh3-seqtk-771d60b/seqtk seq -A - > results/4961-campy_consensus.fasta
# samtools mpileup -u 2638-N12-STEC_V_High_vs_mlst.sorted.bam -f mlstdb/ecoli/baits.fasta |   bcftools view -cg - | vcfutils.pl vcf2fq -d 1| bin/lh3-seqtk-771d60b/seqtk seq -A - > 2638-N12-STEC_V_High_vs_mlst.consensus.fasta

#/mnt/fast/software/blast-2.2.22/bin/blastall -p blastn -i  2638-N12-STEC_V_High_vs_mlst.consensus.fasta -d mlstdb/ecoli/alleles.fas > 2638-N12-STEC_V_High_vs_mlst.consensus.blast.txt

#find . -name "*relreport-176.txt"  | xargs ../bin/metaphlan-1.7.6/merge_metaphlan_tables.py  > merged_table_176.txt

cmdline.run (options)
