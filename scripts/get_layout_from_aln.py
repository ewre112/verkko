#!/usr/bin/env python

import sys
import random

MAX_END_CLIP=50

ref_file = sys.argv[1]
paf_file = sys.argv[2]

def add_lines(read_name, current_lines, lines_per_contig):
	if len(current_lines) == 0: return
	max_mapq_lines = []
	max_mapq = 0
	for l in current_lines:
		parts = l.split('\t')
		# we skip alignments if the full read doesn't align
		# however, in the case of edge effects, we can have reads that extend past the consensus we aligned to and that's OK so allow if the alignment break is within a few bp of contig start/end
		if float(int(parts[3]) - int(parts[2])) / float(int(parts[1])) < 0.95:
			# read is not close to the end so skip it
			if (int(parts[7]) > MAX_END_CLIP and int(parts[8]) + MAX_END_CLIP < int(parts[6])):
				sys.stderr.write("Skipping alignment '%s' because it doesn't cover enough of the read\n"%(l)); continue
		if int(parts[11]) < max_mapq: continue
		if int(parts[11]) > max_mapq:
			max_mapq_lines = []
			max_mapq = int(parts[11])
		if int(parts[11]) == max_mapq: max_mapq_lines.append(l)
	if len(max_mapq_lines) == 0: return
	assert len(max_mapq_lines) >= 1
	# don't allow ambiguously aligned reads
	if len(max_mapq_lines) >= 2: return
	l = max_mapq_lines[0]
	parts = l.split('\t')
	left_clip = int(parts[2])
	right_clip = int(parts[1]) - int(parts[3])
	start_pos = int(parts[7])
	end_pos = int(parts[8])
	if parts[4] == "-":
		end_pos += left_clip
		start_pos -= right_clip
		(start_pos, end_pos) = (end_pos, start_pos)
	else:
		start_pos -= left_clip
		end_pos += right_clip

	lines_per_contig[parts[5]].append((min(int(start_pos), int(end_pos)), start_pos, end_pos, read_name))

lines_per_contig = {}
contig_len = {}
contig_name = ""
current_len = 0
with open(ref_file) as f:
	for l in f:
		if l[0] == '>':
			if contig_name != "":
				lines_per_contig[contig_name] = []
				contig_len[contig_name] = current_len
			contig_name = l[1:].strip()
			current_len = 0
		else:
			current_len += len(l.strip())
if contig_name != "":
	lines_per_contig[contig_name] = []
	contig_len[contig_name] = current_len

current_name = ""
current_lines = []
with open(paf_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] != current_name:
			add_lines(current_name, current_lines, lines_per_contig)
			current_lines = []
			current_name = parts[0]
		current_lines.append(l.strip())

add_lines(current_name, current_lines, lines_per_contig)
current_lines = []
current_name = parts[0]

for contig in lines_per_contig:
	if len(lines_per_contig[contig]) == 0: continue
	lines_per_contig[contig].sort(key=lambda x: x[0])
	print("tig\t" + contig)
	print("len\t" + str(contig_len[contig]))
	print("rds\t" + str(len(lines_per_contig[contig])))
	for line in lines_per_contig[contig]:
		print(line[3] + "\t" + str(line[1]) + "\t" + str(line[2]))
	print("end")
