import scipy.spatial.distance
import numpy
import scipy.cluster.hierarchy
import sys
from Bio import SeqIO

titles = []
def read_contig_stats(fn, method):
	global titles 
        stats = {}
        fh = open(fn)
        cols = fh.readline().rstrip().split("\t")
	titles = cols[1:]

        for ln in fh:
                cols = ln.rstrip().split("\t")
		if method == 'log':
                	stats[cols[0]] = [numpy.log(float(x)+1) for x in cols[1:]]
		elif method == 'count':
			stats[cols[0]] = [float(x) for x in cols[1:]]
		elif method == 'bool':
			stats[cols[0]] = [bool(float(x)) for x in cols[1:]]

#               print stats[cols[0]]
        return stats

distance_method = scipy.spatial.distance.braycurtis
stats_method = 'count'
#stats_method = 'bool'

stats = read_contig_stats(sys.argv[1], method=stats_method)

to_analyse = [rec.id for rec in SeqIO.parse(open(sys.argv[2]), "fasta")]
median_stats = [[] for n in xrange(0, len(stats[to_analyse[0]]))]
for s_n, k in enumerate(to_analyse):
	for n, val in enumerate(stats[k]):
		median_stats[n].append(val)

to_analyse = [numpy.median(n) for n in median_stats]
print ["%.2f" % (t,) for t in to_analyse]

for n, t in enumerate(titles):
	print t, to_analyse[n]

to_test = [rec for rec in SeqIO.parse(open(sys.argv[3]), "fasta")]
for test in to_test:
	t = test.id
	if distance_method(to_analyse, stats[t]) < 0.35:
		SeqIO.write([test], sys.stdout, "fasta")

