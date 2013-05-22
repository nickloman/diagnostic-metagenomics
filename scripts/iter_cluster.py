import pysam
from Bio import SeqIO
from collections import defaultdict
import sys
import itertools
import scipy.spatial.distance
import numpy

#contigs = set([c.rstrip() for c in open("weirdcontigs.txt")])

contigs_to_find = set([rec.id for rec in SeqIO.parse(open(sys.argv[3]), "fasta")])

must_be_in_set = set([rec.id for rec in SeqIO.parse(open(sys.argv[4]), "fasta")])

def read_contig_stats(fn):
	stats = {}
	fh = open(fn)
	fh.readline()
	for ln in fh:
		cols = ln.rstrip().split("\t")
		stats[cols[0]] = [float(x) for x in cols[1:]]
#		print stats[cols[0]]
	return stats

#contig_stats = read_contig_stats(sys.argv[4])

contigs = set([rec.id for rec in SeqIO.parse(open(sys.argv[1]), "fasta")])

f = pysam.Samfile(sys.argv[2], 'rb')
recs = [rec for rec in f if rec.mate_is_unmapped == 0 and rec.alen == 150]

rechash = defaultdict(list)
for rec in recs:
	rechash[rec.qname].append(rec)

recs = [reclist for reclist in rechash.values() if len(reclist) == 2]

for iter_number in xrange(0,100):
	set_found = set()
	path_found = defaultdict(int)

	print >>sys.stderr, "Iter: %d, Total to find: %s, Number found: %s, Correct found: %s, New found: %s, False positives: %s" % (iter_number, len(contigs_to_find), len(contigs), len(contigs & contigs_to_find), len(set_found & contigs_to_find), len(contigs - contigs_to_find))

	for rec in recs:
		pairs = list(itertools.combinations(rec, 2))
		for pair1, pair2 in pairs:
        		contig1 = f.getrname(pair1.tid)
        		contig2 = f.getrname(pair2.tid)

#			print contig1, contig2, pair1.pos, pair1.alen, pair1.aend, pair2.pos, pair2.alen, pair2.aend

#			if contig1 == 'contig-21' and contig2 != 'contig-21':
#				print contig1, contig2
#				print pair1
#				print pair2

    	   	 	if contig1 != contig2:
        	        	if (contig1 in contigs and contig2 not in contigs):
					path_found["%s|%s" % (contig1, contig2)] += 1
				elif (contig2 in contigs and contig1 not in contigs):
					path_found["%s|%s" % (contig2, contig1)] += 1

	#print path_found
	for path, value in path_found.iteritems():
		if value >= 10:
			contig = path.split("|")[1]
			orig_contig = path.split("|")[0]
#			if contig in contigs_to_find:
#				print "GOOD: ",
#			else:
#				print "BAD: ",

#			print "%s > %s: %s" % (contig, orig_contig, scipy.spatial.distance.braycurtis(contig_stats[contig], contig_stats[orig_contig]))
#			print "%s > %s: %s" % (contig, orig_contig, scipy.spatial.distance.dice(contig_stats[contig], contig_stats[orig_contig]))
#			if scipy.spatial.distance.braycurtis(contig_stats[contig], contig_stats[orig_contig]) < 0.15:
#			if scipy.spatial.distance.dice(contig_stats[contig], contig_stats[orig_contig]) < 0.05:
			if contig in must_be_in_set:
				set_found.add(contig)

	print >>sys.stderr, "New found: %s" % (len(set_found & contigs_to_find),)
	contigs = contigs | set_found

	if len(set_found & contigs_to_find) == 0:
		break

print "\n".join(contigs)

