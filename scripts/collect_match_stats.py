import sqlite3
import sys
import os
import re

def get_stats(fn):
        conn = sqlite3.connect(fn)
        conn.row_factory=sqlite3.Row
        c = conn.cursor()
	c.execute("SELECT rowid, * FROM runs")
	rows = c.fetchall()
	outfh = sys.stdout
	for run in rows:
		reads = None
		unaligned = None
		percent = None
		print run['Label']
		try:
			for ln in open(run['Label'] + '.stats.txt'):
				m = re.match("(\d+) reads;", ln)
				if m:
					reads = int(m.group(1))
				m = re.search("(\d+) .* aligned 0 times", ln)
				if m:
					unaligned = int(m.group(1))
			c.execute("UPDATE runs SET TotalReads = ?, HumanMapped = ? WHERE rowid = ?", (reads, reads - unaligned, run['rowid']))
		except Exception, e:
			print e
	conn.commit()
#    outfh.close()

get_stats(sys.argv[1])

