#!/usr/bin/perl

#ASV fasta collection
#open(IN,$ARGV[0]) or die "no $ARGV[0]";
#while($in = <IN>) {
#	if ($in =~ />/) {
#		$in =~ s/>//;
#		chomp $in;
#		$name = $in;
#		$key{$in} = 1;
#		print STDERR "key $name found\n";
#	} else {
#		$seq{$name} = $in;
#	}
#}
#close IN;

#print STDERR scalar keys %key, " keys collected\n"; 

#Bowtie output
open(IN,$ARGV[0]) or die "no $ARGV[0]";
$asv = 0;
while($in = <IN>) {
	chomp $in;
	next if ($in !~ /\tchr/);
	@tmp = split ("\t",$in);
	$chrom{$tmp[0]} = $tmp[2];
	$pos{$tmp[0]} = $tmp[3];
	$seq{$tmp[0]} = $tmp[9];
#	print STDERR "$tmp[0] chrom:$tmp[2] pos:$tmp[3] seq:$tmp[9]\n";
	$hseq++;
}
close IN;

print STDERR scalar keys %seq," human sequences found by bowtie\n"; 

#count collection (a specific OTU/ASV table with ASV IDs as row names)
$totASV = -1;
$hASV = 0;
$totseq = 0;
$hcount = 0;
$totCount = 0;
$totCounts = 0;
$totASV = 0;
$out = '';

#OTU TABLE
open(IN,$ARGV[1]) or die "no $ARGV[1]";
@in = <IN>;
close IN;
shift @in; #first line are headers, removed.
$totASV = scalar @in;
foreach $in (@in) {
	chomp $in;
	@tmp = split (/[\t;]/,$in); #elements are split
	$ASVname = shift @tmp; #row name (ASV name) is removed from array and stored
	splice @tmp, -7; #last 7 are taxonomy, rmeoved
	#<STDIN>;
	$nsamp = scalar @tmp -1; #number of samples (1st is the ASV name
	#<STDIN>;
	$ASVcount = eval(join "+", @tmp); #grand reads total of the current row
	$ASVcounts += $ASVcount; #grand reads total of the current row
	if ($seq{$ASVname}) { #this menas a match with bowtie
		#print STDERR "key $tmp[0] found\n";
		$hASV++; #counter for humber of human ASV
		$hASVcounts += $ASVcount; #counter for abundance of human ASV
		$nsamp = 0;
		$hsamp = 0;
		foreach my $v (@tmp) {
			$hsamp++ if ($v > 0);
			$nsamp++;
		}
		$out .= "$ASVname\t$ASVcount\t$hsamp\t$nsamp\t".$chrom{$ASVname}."\t".length($seq{$ASVname})."\t".$seq{$ASVname}."\n" if ($hASVcounts>0);
		$chrcount{$chrom{$ASVname}} += $ASVcount;
#		print STDERR "$ASVname\t$hASVcounts\t$hsamp\t$nsamp\t".$chrom{$tmp[0]}."\t".length($seq{$tmp[0]})."\t".$seq{$tmp[0]}."\n" if ($tot>0);
		#print STDERR "H: $ASVname in $hsamp / $nsamp samples, count $ASVcount\n";
	} else {
#		print STDERR "$ASVname: $nsamp samples, $ASVcount total counts, NON human\n";
		#$totcount += $ASVcount; #grand total of all matched hASV
		#push(@no,$tmp[0]);
	}
	#<STDIN>;
}

close IN;

print STDERR "#ASV\thASV\thcount\ttotcount\n#$totASV\t$hASV\t$hASVcounts\t$ASVcounts\n";
print "#ASV\thASV\thcount\ttotcount\n#$totASV\t$hASV\t$hASVcounts\t$ASVcounts\n";
foreach $k (sort {$chrcount{$b} <=> $chrcount{$a} } keys %chrcount) {
	print "#$k\t$chrcount{$k}\n";
}
print "ID\tcount\tnsamples\tchrom\tseqlen\tseq\n$out\t";;
#print STDERR "ASV\thASV\thcount\ttotCount\n$totASV\t$hASV\t$tot\t$totCount\n";

