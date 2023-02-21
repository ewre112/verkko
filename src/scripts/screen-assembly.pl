#!/usr/bin/env perl

use strict;

my $assembly   = "";
my %contigLength;
my $graph      = "";
my $graphmap   = "";
my %gtoc;
my $hificov    = "";
my %hifiCoverage;
my $minLength  = 100000;
my %contaminantSeq;

my $threads    = 4;
my $mashmap    = "mashmap";

my $output;

while (scalar(@ARGV) > 0) {
    my $opt = shift @ARGV;

    if    ($opt eq "--assembly") {
        $assembly = shift @ARGV;
    }
    elsif ($opt eq "--graph") {
        $graph = shift @ARGV;
    }
    elsif ($opt eq "--graphmap") {
        $graphmap = shift @ARGV;
    }
    elsif ($opt eq "--hifi-coverage") {
        $hificov = shift @ARGV;
    }
    elsif ($opt eq "--minlength") {
        $minLength = int(shift @ARGV);
    }
    elsif ($opt eq "--contaminant") {
        my $prefix   = $ARGV[0];
        my $sequence = $ARGV[1];

        while (($prefix !~ m/^--/) && (-e $sequence)) {
            $contaminantSeq{$prefix} = $sequence;

            shift @ARGV;
            shift @ARGV;

            $prefix   = $ARGV[0];
            $sequence = $ARGV[1];
        }
    }
    elsif ($opt eq "--threads") {
        $threads = int(shift @ARGV);
    }
    elsif ($opt eq "--output") {
        $output = shift @ARGV;
    }
    elsif ($opt eq "--mashmap") {
        $mashmap = shift @ARGV;
    }
    else {
        die "Unknown option '$opt'.\n";
    }
}

if (($assembly eq "") ||
    ($graph    eq "") ||
    ($output   eq "")) {
    print "usage: $0 --assembly X.fasta --graph X.gfa ...\n";
    print "  --assembly X.fasta           sequences to screen for crud\n";
    print "  --graph X.gfa                graph to decide if a short contig is disconnected\n";
    print "  --graphmap X.scfmap          layoutContigs scfmap to rename graph nodes to contig names\n";
    print "  --hifi-coverage X.csv        contig read coverage\n";
    print "  --minlength L                contigs shorter than L bp are deemed 'short' (default: 100000)\n";
    print "  --contaminant N F [N F ...]  label N and contaminant fasta F to screen; example:\n";
    print "                                  --contaminant ebv  ebv.fasta.gz \\\n";
    print "                                                mito mito.fasta   \\\n";
    print "                                                rdna r.fasta.gz\n";
    print "  --threads T                  number of mashmap compute threads to use\n";
    print "  --output X                   write outputs to prefix X\n";
    print "                                  X.disconnected.fasta\n";
    print "                                  X.ebv.exemplar.fasta\n";
    print "                                  X.ebv.fasta\n";
    print "                                  X.ebv.mashmap.err\n";
    print "                                  X.ebv.mashmap.out\n";
    print "                                  X.fasta\n";
    print "                                  X.mito.exemplar.fasta\n";
    print "                                  X.mito.fasta\n";
    print "                                  X.mito.mashmap.err\n";
    print "                                  X.mito.mashmap.out\n";
    print "                                  X.rdna.exemplar.fasta\n";
    print "                                  X.rdna.fasta\n";
    print "                                  X.rdna.mashmap.err\n";
    print "                                  X.rdna.mashmap.out\n";
    exit(1);
}





#  Read sequence lengths into global %contigLength.
#
sub loadSequenceLengths ($) {
    my $a = shift @_;
    my $n = 0;
    my $t = 0;

    print STDERR "Loading sequence lengths from '$a'.\n";

    open(F, "< $a") or die "Failed to open '$a' for reading: $!\n";
    while (!eof(F)) {
        my $h = <F>;
        my $s = <F>;

        $h =~ s/^\s+//;   $s =~ s/^\s+//;
        $h =~ s/\s+$//;   $s =~ s/\s+$//;

        if ($h =~ m/^>(\S+)\s*/) {
            $contigLength{$1}  = length($s);
            $n                += 1;
            $t                += length($s);
        }
        else {
            die "Failed to find ident line in input line '", substr($h, 0, 40), "'.\n";
        }
    }
    close(F);

    print STDERR "  Found $n sequences with total length $t bp.\n";
    print STDERR "\n";
}



#  Read hifi coverage for the contigs.
#
sub loadHifiCoverage ($) {
    my $c = shift @_;

    print STDERR "Loading contig coverage from '$c'.\n";

    open(F, "< $c") or die "Failed to open '$c' for reading: $!.\n";
    while (<F>) {
        my ($ident, $cov) = split '\s+', $_;

        if (exists($gtoc{$ident})) {
            #print STDERR "$ident -> $gtoc{$ident} -> $cov\n";

            $hifiCoverage{$gtoc{$ident}} = $cov;
        }
        else {
        }
    }
    close(F);

    print STDERR "  Found ", scalar(keys %hifiCoverage), " contig coverage values.\n";
    print STDERR "\n";
}



#  Load a map from graph-contig-names to assembly-contig-names.
#
sub loadGtoC ($) {
    my $m = shift @_;

    print STDERR "Loading graph to contig name map from '$m'.\n";

    open(M, "< $m") or die "Failed to open graph-to-contig name map '$m' for reading: $!\n";
    while (<M>) {
        s/^\s+//;
        s/\s+$//;

        if (m/^path\s+(\S+)\s+(\S+_from_|\S+_unused_){0,1}(\S+)$/) {
            my $cn = $1;
            my $gn = $3;
            $gtoc{$gn} = $cn;
        } elsif (m/path/) {
            die "Failed to parse path '$_'\n";
        }
    }
    close(M);

    print STDERR "  Found ", scalar(keys %gtoc), " contig names.\n";
    print STDERR "\n";
}



#  Find disconnected nodes in the graph.  Returns an array of node names.
#
#  cat ../{input.graph} | grep ^L |                 awk '{{ print \$2"\n"\$4 }}' | sort | uniq > ./nodes-with-links
#  cat ../{input.graph} | grep ^S | sed s/LN:i:// | awk '{{ print \$2        }}' | sort | uniq > ./nodes
#
#  diff --side-by-side --suppress-common-lines ./nodes ./nodes-with-links | awk '{{ print \$1 }}' > ./nodes-isolated
#
sub findDisconnected ($$) {
    my $g = shift @_;
    my $minl = shift @_;
    my %nodes;
    my %edges;
    my @discShort;
    my @discLong;

    print STDERR "Finding disconnected nodes in '$g' shorter than $minl.\n";

    open(G, "< $g") or die "Failed to open graph '$g' for reading: $!\n";
    while (<G>) {
        my ($t, $n1, undef, $n2) = split '\t', $_;

        #next   if (!exists($gtoc{$n1}));   #  The graph seems to have stuff not in the output?
        #next   if (!exists($gtoc{$n2}));

        #  Not sure this is correct.  On a rukki run, the warnings below
        #  trigger for some nodes.  I suspect these are the things rukki is
        #  building scaffolds with.

        if    ($t eq "S") {
            my $c1 = $gtoc{$n1};

            #print STDERR "WARNING1: graph name '$n1' not found.\n"   if (!defined($c1));

            $nodes{$c1}++;
        }
        elsif ($t eq "L") {
            my $c1 = $gtoc{$n1};
            my $c2 = $gtoc{$n2};

            #print STDERR "WARNING2: graph name '$n1' not found.\n"   if (!defined($c1));
            #print STDERR "WARNING2: graph name '$n2' not found.\n"   if (!defined($c2));

            $edges{$c1}++;
            $edges{$c2}++;
        }
    }
    close(G);

    foreach my $c (keys %nodes) {
        next   if (exists($edges{$c}));

        if ($contigLength{$c} <  $minl) {
            push @discShort, $c;
        } else {
            push @discLong, $c;
        }
    }

    printf STDERR " Scanned %7d nodes.\n", scalar(keys %nodes);
    printf STDERR " Scanned %7d edges.\n", scalar(keys %edges);
    printf STDERR " Found %5d short disconnected nodes.\n", scalar(@discShort);
    printf STDERR " Found %5d long  disconnected nodes.\n", scalar(@discLong);
    printf STDERR "\n";

    return(@discShort);
}



sub mashMap ($$$$) {
    my $a      = shift @_;
    my $output = shift @_;
    my $cp     = shift @_;
    my $cs     = shift @_;

    my $map;

    if (-e "$output.$cp.mashmap.out") {
        print STDERR "Using pre-computed mashmap results for '$cp' ('$cs').\n";
    }
    else {
        print STDERR "Running mashmap against '$cp' ('$cs').\n";

        $map  = "$mashmap";
        $map .= " --ref '$a'";
        $map .= " --query '$cs'";
        $map .= " --perc_identity 95";
        $map .= " --segLength 10000";
        $map .= " --filter_mode none";
        $map .= " --threads $threads";
        $map .= " --output '$output.$cp.mashmap.out'";
        $map .= " > $output.$cp.mashmap.err 2>&1";

        system($map);
    }
}



#
#  Filters mashmap output to remove low-quality and spurious matches.
#  Remembers the best match.
#
#    pick node that contains most reference
#    and has highest depth
#      pick the five/ten highest coverage contigs
#      then pick the highest depth
#
#
#  MashMap output columns:
#     1 query
#     2 len
#     3 bgn
#     4 end
#     5 strand
#     6 ref
#     7 len
#     8 bgn
#     9 end
#    10 ident
#
#  awk '{{if (\$NF > 99 && ((\$9-\$8)/\$7 > 0.50 || (\$9-\$8)/\$7 > 0.25 && \$7 < 50000)) print \$6}}' | sort | uniq > contaminant.list
#
sub filterMash ($) {
    my  $cp = shift @_;

    #
    #  Pass 1: find the highest reference contig coverage.  The 'reference contig' is
    #  the first sequence...called 'query' above.
    #

    my $qcovBest = 0;   #  Best qcov encountered.
    my $qcovMax  = 0;   #  Number of times we've seen > 99% qcov.

    open(F, "< $output.$cp.mashmap.out") or die "Failed to open mashmap output '$output.$cp.mashmap.out' for reading: $!\n";
    while (<F>) {
        my ($q, $qlen, $qbgn, $qend, $s, $r, $rlen, $rbgn, $rend, $ident) = split '\s+', $_;

        my $qcov = int(0.5 + 100.0 * ($qend - $qbgn) / $qlen);
        my $rcov = int(0.5 + 100.0 * ($rend - $rbgn) / $rlen);

        my $gi   = ($ident >= 99.0);
        my $g1   = ($rcov  >= 50);
        my $g2   = ($rcov  >= 25) && ($rlen < 50000);

        if (($gi) && ($g1 || $g2)) {                      #  If ident is good and the contig
            $qcovBest = $qcov   if ($qcovBest <  $qcov);  #  is primarily the contaminant,
            $qcovMax += 1       if ($qcov     >= 100);    #  save the best contaminant coverage.
        }
    }

    #  Pick a minimum acceptable qcov for filtering.  If there is anything with full coverage,
    #  only use full coverage, otherwise, take the best found and anything close to it.
    my $qcovMin = ($qcovMax > 0) ? (100) : ($qcovMax - 5);

    #print "qcovMin=$qcovMin  qcovBest=$qcovBest  qcovMax=$qcovMax\n";

    #
    #  Pass 2: find the exemplar.
    #

    my $nHits;
    my $nSaveI;
    my $nSave1;
    my $nSave2;

    my %hits;

    my $bestIdent = "";
    my $bestDepth = 0;

    open(F, "< $output.$cp.mashmap.out") or die "Failed to open mashmap output '$output.$cp.mashmap.out' for reading: $!\n";
    while (<F>) {
        my ($q, $qlen, $qbgn, $qend, $s, $r, $rlen, $rbgn, $rend, $ident) = split '\s+', $_;

        my $qcov = int(0.5 + 100.0 * ($qend - $qbgn) / $qlen);
        my $rcov = int(0.5 + 100.0 * ($rend - $rbgn) / $rlen);

        my $gi   = ($ident >= 99.0);
        my $g1   = ($rcov  >= 50);
        my $g2   = ($rcov  >= 25) && ($rlen < 50000);
        my $g3   = ($qcov  >= $qcovMin);

        $nHits++;
        $nSaveI += 1   if ($gi);
        $nSave1 += 1   if ($g1) && (!$gi);
        $nSave2 += 1   if ($g2) && (!$gi);

        if (($gi) && ($g1 || $g2)) {            #  If ident is good and the contig is primarily
            $hits{$r}++;                        #  contaminant, save the id for filtering.
 
            #print STDERR "bestDepth='$bestDepth'  r='$r'  hc='$hifiCoverage{$r}'\n";

            if (($g3) && ($bestDepth <= $hifiCoverage{$r})) {
                $bestIdent = $r;
                $bestDepth = $hifiCoverage{$r};
            }
        }
    }
    close(F);

    printf "  Found %5d       hits.\n", $nHits;
    printf "  Found %5d       hits with good identity.\n", $nSaveI;
    printf "  Found %5d       hits with good coverage but bad identity.\n", $nSave1;
    printf "  Found %5d short hits with good coverage but bad identity.\n", $nSave2;
    printf "  Found %5d       contaminant contigs.\n", scalar(keys %hits);
    printf "  Best hit has %6.2fx read coverage. ('%s').\n", $bestDepth, $bestIdent;
    printf "\n";

    return($bestIdent, keys %hits);
}



#  Extract '@filter' sequences from $a.
#
#  If $cp is defined, the sequences in @filter are copied to the output file.
#  If it is undef, the sequences NOT in @filter are copied to the output file.
#
#  If $exemplar is defined, that specific sequence is saved.
#
sub filterSequences ($$$$@) {
    my $a           = shift @_;
    my $output      = shift @_;
    my $cp          = shift @_;
    my $exemplar    = shift @_;
    my @filter      =       @_;
    my %filter;

    #  Make lookups of sequences we want to filter easier.
    open(O, "> $output.ids") or die "Failed to open '$output.ids' for output: $!\n";
    foreach my $f (@filter) {
        print O "$f\n";
        $filter{$f}++;
    }
    close(O);

    #  Open some output files.
    if (defined($cp)) {
        print STDERR "Filtering '$cp' from '$a'.\n";
        open(O, "> $output.fasta")          or die "Failed to open '$output.fasta' for output: $!\n";
        open(E, "> $output.exemplar.fasta") or die "Failed to open '$output.exemplar.fasta' for output: $!\n";
    }
    else {
        print STDERR "Extracting sequences from '$a'.\n";
        open(O, "> $output.fasta")          or die "Failed to open '$output.fasta' for output: $!\n";
    }

    #  Scan the input, copying desired sequences to the output.
    open(F, "< $a") or die "Failed to open '$a' for reading: $!\n";
    while (!eof(F)) {
        my $h = <F>;
        my $s = <F>;

        $h =~ s/^\s+//;   $s =~ s/^\s+//;
        $h =~ s/\s+$//;   $s =~ s/\s+$//;

        if ($h =~ m/^>(\S+)\s*/) {
            my $n = $1;

            if (defined($cp)) {
                print O "$h\n$s\n"  if (exists($filter{$n}));
                print E "$h\n$s\n"  if ($n eq $exemplar);
            }
            else {
                print O "$h\n$s\n"  if (!exists($filter{$n}));
            }
        }
        else {
            die "Failed to find ident line in input line '", substr($h, 0, 40), "'.\n";
        }
    }
    close(F);

    close(E)   if (defined($cp));
    close(O);

    unlink "$output.exemplar.fasta" if (!defined($exemplar));

    print STDERR "\n";
}






#
#  Main
#

loadSequenceLengths($assembly);
loadGtoC($graphmap);

my @disconnected  = findDisconnected($graph, $minLength);
filterSequences($assembly, "$output.disconnected", "disconnected", undef, @disconnected);

my @crud = @disconnected;

loadHifiCoverage($hificov);

foreach my $cp (keys %contaminantSeq) {
    mashMap($assembly, $output, $cp, $contaminantSeq{$cp});

    my ($exemplar, @contaminants) = filterMash($cp);

    filterSequences($assembly, "$output.$cp", $cp, $exemplar, @contaminants);

    push @crud, @contaminants;
}


filterSequences($assembly, $output, undef, undef, @crud);

exit(0);
