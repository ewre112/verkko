MBG -i hifi.fa -o hifi-resolved.gfa -k 1001 -w 100 -a 1 -u 2 --error-masking=collapse-msat --output-sequence-paths paths.gaf -r 15000

scripts/insert_aln_gaps.py hifi-resolved.gfa 2 50 gaps-hifi-1.gaf gapone < paths.gaf > gapped-once-hifi-resolved.gfa
scripts/insert_aln_gaps.py gapped-once-hifi-resolved.gfa 2 300 gaps-hifi-2.gaf gaptwo < paths.gaf > gapped-twice-hifi-resolved.gfa
scripts/insert_aln_gaps.py gapped-twice-hifi-resolved.gfa 1 5 gaps-hifi-3.gaf gapthree < paths.gaf > gapped-hifi-resolved.gfa
cut -f 6 < paths.gaf | scripts/unroll_tip_loops.py gapped-hifi-resolved.gfa 5 > unrolled-hifi-resolved.gfa
scripts/get_unroll_mapping.py gapped-hifi-resolved.gfa unrolled-hifi-resolved.gfa > unroll_mapping_1.txt
scripts/unitigify.py "utig1-" unitig-mapping-1.txt < unrolled-hifi-resolved.gfa > unitig-unrolled-hifi-resolved.gfa

# only for evaluating hifi-only CHM13 haploid assemblies
# scripts/untip_relative.py 15000 15000 0.1 0.1 < unitig-unrolled-hifi-resolved.gfa > hifi-resolved-graph-tip.gfa
# scripts/pop_bubbles_keep_longest.py 10 < hifi-resolved-graph-tip.gfa > popped-hifi-resolved-graph-tip.gfa
# scripts/unitigify.py < popped-hifi-resolved-graph-tip.gfa > unitig-popped-hifi-resolved-graph-tip.gfa

GraphAligner -t 32 -g unitig-unrolled-hifi-resolved.gfa -f ont.fa -a alns-ont.gaf --seeds-mxm-length 30 --seeds-mem-count 10000 -b 15 --multimap-score-fraction 0.99 --precise-clipping 0.7 --min-alignment-score 5000 --hpc-collapse-reads --discard-cigar 1> stdout_ga_ont.txt 2> stderr_ga_ont.txt

awk -F '\t' '{if ($4-$3 >= $2*0.8 && $12 >= 20) print;}' < alns-ont.gaf > alns-ont-filter.gaf
scripts/trim_dbg_alignment.py unitig-unrolled-hifi-resolved.gfa 1500 < alns-ont-filter.gaf > alns-ont-filter-trim.gaf
scripts/calculate_coverage.py unitig-unrolled-hifi-resolved.gfa < alns-ont-filter-trim.gaf > nodecovs-ont.csv

cut -f 6 < alns-ont-filter-trim.gaf > paths.txt
awk -F '\t' '{if ($12 >= 20) print;}' < alns-ont.gaf > alns-ont-mapqfilter.gaf
scripts/insert_aln_gaps.py unitig-unrolled-hifi-resolved.gfa 3 50 gaps-ont.gaf gapont < alns-ont-mapqfilter.gaf > gapped-unitig-unrolled-hifi-resolved.gfa
awk '{if ($2 >= 100000) {sum += $2*$3; count += $2;}}END{print sum/count;}' < nodecovs-ont.csv
scripts/estimate_unique_local.py gapped-unitig-unrolled-hifi-resolved.gfa alns-ont-filter-trim.gaf 100000 30 0.8 > unique_nodes_ont_coverage.txt
# scripts/translate_uniques.py normal-hifi_connected_twice.gfa < unique_nodes_hifi.txt > translated_uniques.txt
# scripts/translate_nodes_by_seq.py normal-hifi_connected_twice.gfa unitig-unrolled-hifi-resolved.gfa < translated_uniques.txt > unique_nodes_ont_translated.txt
# cat unique_nodes_ont_coverage.txt unique_nodes_ont_translated.txt | sort | uniq > unique_nodes_ont.txt
scripts/fix_diploid_unique_nodes.py unique_nodes_ont_coverage.txt nodecovs-ont.csv gapped-unitig-unrolled-hifi-resolved.gfa > unique_nodes_diploidfix.txt
cp unique_nodes_diploidfix.txt unique_nodes_ont.txt

scripts/find_bridges.py unique_nodes_ont.txt < paths.txt > bridges.txt
grep -v '(' < bridges.txt | grep -vP '^$' | scripts/remove_wrong_connections_2.py forbidden_wrong_connections.txt | sort > bridging_seq_all.txt
scripts/pick_majority_bridge.py forbidden_minority_bridges.txt < bridging_seq_all.txt > bridging_seq_picked_all.txt
scripts/remove_crosslink_paths.py unique_nodes_ont.txt bridging_seq_picked_all.txt bridges.txt > bridges_fixcrosslink.txt 2> forbidden_crosslinks.txt
scripts/fix_diploid_paths.py unique_nodes_ont.txt gapped-unitig-unrolled-hifi-resolved.gfa bridges_fixcrosslink.txt bridges.txt 3 > bridging_seq_diploidfix_all.txt
cp bridging_seq_diploidfix_all.txt bridging_seq_picked.txt
# forbidden_wrong_connections.txt deliberately not included here so that if that causes a gap, the tangle is forbidden
cat forbidden_crosslinks.txt forbidden_minority_bridges.txt > bridging_seq_forbidden.txt
scripts/forbid_unbridged_tangles.py unique_nodes_ont.txt gapped-unitig-unrolled-hifi-resolved.gfa bridging_seq_forbidden.txt bridging_seq_picked.txt paths.txt nodecovs-ont.csv 30 > forbidden_ends.txt
scripts/connect_uniques.py gapped-unitig-unrolled-hifi-resolved.gfa forbidden_ends.txt bridging_seq_picked.txt > connected.gfa

scripts/merge_unresolved_dbg_nodes.py < connected.gfa > normal-connected.gfa
scripts/get_bridge_mapping.py normal-connected.gfa gapped-unitig-unrolled-hifi-resolved.gfa > bridge_mapping.txt
scripts/add_fake_alignments.py unitig-unrolled-hifi-resolved.gfa normal-connected.gfa alns-ont-filter-trim.gaf nodecovs-ont.csv fake-ont-alns.gaf fake-ont-nodecovs-once.csv 10
scripts/add_fake_bridging_paths.py forbidden_ends.txt bridging_seq_picked.txt fake-ont-nodecovs-once.csv fake-ont-nodecovs.csv 10 >> fake-ont-alns.gaf
/usr/bin/time -v scripts/resolve_triplets_kmerify.py normal-connected.gfa fake-ont-paths.txt fake-ont-nodecovs.csv resolve-mapping.txt 100000 3 5 3 2 < fake-ont-alns.gaf > ont-resolved-graph.gfa 2> stderr_ont_resolved_graph.txt

scripts/unroll_tip_loops.py ont-resolved-graph.gfa 3 < fake-ont-paths.txt > unrolled-ont-resolved.gfa
scripts/get_unroll_mapping.py ont-resolved-graph.gfa unrolled-ont-resolved.gfa > unroll_mapping_2.txt
scripts/unitigify.py "utig2-" unitig-mapping-2.txt < unrolled-ont-resolved.gfa > unitig-unrolled-ont-resolved.gfa

# consensus
grep -P '^S' < unitig-unrolled-ont-resolved.gfa | awk '{print ">" $2; print $3;}' > contigs_rle.fa
winnowmap -x map-pb -t 32 contigs_rle.fa hifi.fa > alns.paf
cat gaps-ont.gaf | cut -f 1 | sort | uniq > used_ont.txt
scripts/pick_reads_stdin.py used_ont.txt < ont.fa > ont_gap_subset.fa
scripts/rle.py < ont_gap_subset.fa > ont_gap_subset_rle.fa
winnowmap -x map-ont -t 32 contigs_rle.fa ont_gap_subset_rle.fa >> alns.paf
scripts/get_layout_from_aln.py contigs_rle.fa alns.paf read_names.txt hifi.fa ont_gap_subset.fa > layout.txt

# layout without alignment
cat *mapping* > combined-nodemap.txt
cat *.gfa | grep -P '^L' > combined-edges.gfa
cat gaps-*.gaf paths.gaf > combined-alignments.gaf
# remove * so that noseq-..gfa don't mess up the node lengths
grep -P '^S' *.gfa | grep -v '\*' | awk '{print $2 "\t" length($3);}' > nodelens.txt
grep -P '^S' unitig-unrolled-ont-resolved.gfa | awk '{print $2 "\t" ">" $2;}' > consensus_paths.txt
scripts/get_layout_from_mbg.py combined-nodemap.txt combined-edges.gfa combined-alignments.gaf consensus_paths.txt nodelens.txt > layout.txt 2> unitig_to_mbg_list.txt

# just for debug info
scripts/check_layout_gaps.py < layout.txt > gaps.txt

# only for evaluating CHM13 haploid assemblies
# scripts/untip_relative.py 30000 30000 0.1 0.1 < unitig-unrolled-ont-resolved.gfa > connected-tip.gfa
# scripts/unitigify.py "utig3-" unitig-mapping-3.txt < connected-tip.gfa > unitig-normal-connected-tip.gfa
# scripts/pop_bubbles_keep_longest.py 10 < unitig-normal-connected-tip.gfa > popped-unitig-normal-connected-tip.gfa
# scripts/unitigify.py "utig4-" unitig-mapping-4.txt < popped-unitig-normal-connected-tip.gfa > unitig-popped-unitig-normal-connected-tip.gfa
