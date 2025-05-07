#!/usr/bin/perl

#ASV collection
open(IN,$ARGV[0]) or die "no";
while($in = <IN>) {
	if ($in =~ />/) {
		$in =~ s/>//;
		chomp $in;
		$key{$in} = 1;
	} else {
		$seq{$in} = $in;
	}
}
close IN;

print STDERR scalar keys %key, " keys collected\n"; 

#count collection
open(IN,$ARGV[1]) or die "no";
while($in = <IN>) {
	@tmp = split (/\t/,$in);
	if ($key{$tmp[0]}) {
		#print STDERR "key $tmp[0] found";
		$tot = 0;
		foreach $p (1..$#tmp) {
			#print STDERR $p,"+";
			next if $p =~ /a-zA-Z/;
			$tot += $tmp[$p];
		}
		print $tmp[0],"\t",$tot,"\t",$seq{$in},"\n";
	} else {
		push(@no,$tmp[0]);
	}
}
close IN;

#print "MISSING: ",join "\n",@no,"\n";

