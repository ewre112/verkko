#!/usr/bin/env python

import sys

# input partial layouts from argv[1:]
# output full layout to stdout

contig_lines = {}
contig_len = {}

current_contig = ""
for filename in sys.argv[1:]:
	with open(filename) as f:
		for l in f:
			if l.strip() == "end":
				current_contig = ""
				continue
			parts = l.strip().split('\t')
			if parts[0] == "tig":
				current_contig = parts[1]
			elif parts[0] == "len":
				length = int(parts[1])
				assert current_contig not in contig_len or contig_len[current_contig] == length
				contig_len[current_contig] = length
			elif parts[0] == "rds":
				continue
			else:
				assert current_contig != ""
				contig_lines[current_contig].append((parts[0], int(parts[1]), int(parts[2])))

for contig in contig_lines:
	if len(contig_lines[contig]) == 0: continue
	contig_start = contig_lines[contig][0][1]
	contig_end = contig_lines[contig][0][1]
	for line in contig_lines[contig]:
		contig_start = min(contig_start, line[1])
		contig_start = min(contig_start, line[2])
		contig_end = max(contig_end, line[1])
		contig_end = max(contig_end, line[2])
	print("tig\t" + contig)
	print("len\t" + str(contig_end - contig_start))
	print("rds\t" + str(len(contig_lines[contig])))
	contig_lines[contig].sort(key=lambda x: min(x[1], x[2]))
	for line in contig_lines[contig]:
		assert line[1] >= contig_start
		assert line[2] >= contig_start
		assert line[1] - contig_start <= contig_end
		assert line[2] - contig_start <= contig_end
		print(line[0] + "\t" + str(line[1] - contig_start) + "\t" + str(line[2] - contig_start))
	print("end")
