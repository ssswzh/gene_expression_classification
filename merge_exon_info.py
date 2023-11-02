#!/usr/bin/env python3

import sys
from collections import OrderedDict

if len(sys.argv) == 1 :
	sys.exit("""Usage:
	python {} gencode.v38lift37.annotation.extract_exon.merge.anno.bed > gencode.v38lift37.annotation.extract_exon.merge.anno.final.bed
	""".format(sys.argv[0]) )

bed = open(sys.argv[1])
result = OrderedDict()
for line in bed:
	mchr, mstart, mend, ochr, ostart, oend, gene, exon, strand, length, geneid, transid, exonid, overlaplen = line.strip().split("\t")
	key = ':'.join([mchr, mstart, mend])
	if key in result:
		oldgene, oldstrand, oldlength, oldgeneid, oldtransid, oldexonid = result[key]
		transid = oldtransid + "|" + transid
		exonid = oldexonid + "|" + exonid
		result[key] = [oldgene, oldstrand, oldlength, oldgeneid, transid, exonid]
	else:
		newlength = str(int(mend) - int(mstart))
		result[key] = [gene, strand, newlength, geneid, transid, exonid]

print('\t'.join(['#chrom', 'start', 'end', 'gene', 'strand', 'length', 'gene_id', 'transcript_id', 'exon_id']))
for key in result:
	info = list(key.split(':')) + list(result[key])
	print('\t'.join(info))
