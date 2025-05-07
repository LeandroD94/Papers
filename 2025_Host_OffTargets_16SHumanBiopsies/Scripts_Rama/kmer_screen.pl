#!/usr/bin/perl

if (!$ARGV[0]) {
	print STDERR "USAGE: kmer_screen.pl reads.fasta kmer_size best_kmer_number\n\n";
	exit;
}
open(IN,$ARGV[0]);
while ($s = <IN> ) {
	chomp $s;
	if ($s =~ />/) {
		$name = $s;
	} else {
		$seq{$name} .= $s;
	}
}
close IN;
print STDERR scalar keys %seq," loaded\n";
$ksize = $ARGV[1];
foreach $n (keys %seq) {
	$cnt = 0;
#	print STDERR "Scanning $n\n";
	while($seq{$n} =~ /(.{$ksize})/g) {
		$counts{$1}++;
		pos($seq{$n}) -= $k-1;
		#print STDERR $1,"\n";
		$cnt++;
	}
}
print "kmer\ttotal_mapped\tTm\treads_mappe\ttarget_reads\n";
$max = $ARGV[2];
$cnt = 0;
foreach $kmer (sort {$counts{$b} <=> $counts{$a}} keys %counts) {
	$cnt++;
	last if $cnt > $max;
	foreach $k (keys %seq) {
		if ($seq{$k} =~ /$kmer/) {
			push(@{$targets{$kmer}},$k);	
		}
	}
	@tm = split("\t",`FAStm.pl -$kmer`);
	$tm[1] =~ s/\n//g;
	print $kmer,"\t",$counts{$kmer},"\t",$tm[1],"\t",scalar keys @{$targets{$kmer}},"\t",(join "|",@{$targets{$kmer}}),"\n";
}

