#!/usr/bin/perl
if ($ARGV[0] eq "-h") {
	print "\nUSAGE: FAStm.pl [[-h]123] [-primer|fastafile[.gz]]\n\n";
	print STDERR " -primer: a single primer\n";
	print STDERR " -fastafile: a fasta file containing primer sequences\n";
	print STDERR " (1) 69.3 + 0.41*GC% - 650/L (default)\n (2) 2*AT + 4*GC\n (3) as in primer3, with [2+ ions] = 2.5 mM and [dNTP] = 0.2 nM\n\n";
	exit;
}

$ions2 = 2.5; #divalent ions [mM]
$dntp = 0.2; # primers [nM]

#print STDERR "\n\n --- HIGHLY EXPERIMENTAL: THIS IS TO BE TESTED YET ---\n\n";
#efault mode
$mode = 3;

#if a file is provided it must be in the 1st or last ARGV 
$file = shift @ARGV if (-e $ARGV[0]);
$file = pop @ARGV if (-e $ARGV[$#ARGV]);

#set the Tm calculation mode (see below)
$mode = shift @ARGV if ($ARGV[0] =~ /^\d+$/);
$mode = pop @ARGV if (-d $ARGV[$#ARGV] =~ /^\d+$/);

#what is left can be nothing (STDIN) or preceded by a "-" sign
#the "-" sign marks a deviation form the normal FASsuite operation:
if ($ARGV[0] =~ /-/) {
	$ARGV[0] =~ s/-//;
	#there are primers in a file as a newline separated list
	@tmp = split (/-/,$ARGV[0]);
	foreach $p (@tmp) { 
		$tm = scan($p,$mode);
		print $p,"\t",$tm,"\n";
	}
	exit;
}
#this is the normal FASsuite flow: the file can be a gz file or not, 
#in any case it is redirected to STDIN to allow piping
if($file) {
	open($in,"zcat $file|") if (-e $file && $file =~ /\.gz$/);
	open($in,"$file") if (!-e $file && $file !~ /\.gz$/);
} else {
	$in = *STDIN;	
}
while($l = <$in>) {
	chomp $l;
	if ($l =~ />/) {
		if ($seq) {
			$tm = scan($seq,$mode);
			print $name,"\t",$seq,"\t",$tm,"\n";
		}
		$name = $l;
		$seq = '';
	} else {	
		$seq .= $l;
	}
}

$tm = scan($seq,$mode);
print $name,"\t",$seq,"\t",$tm,"\n";

sub scan {
	my $tm;
	my $seq = shift;
	my $mode = shift;
	$seq =~ s/[\n\r\s]//g;
	$gc = $seq =~ s/([gc])/$1/gi;
	$len = length($seq);
	$gcp = sprintf ("%.2f",$gc/$len*100);
	if ($mode == 1) {
		$tm = 69.3 + (41*$gc/$len)-(650/$len);
		$tm = sprintf("%.2f", $tm);
	}
	#$tm = 81.5 + (16.6*log($Na/(1+0.7*$Na))/log(10))+(41*$gc/$len)-(500/$len); #as seen in the melting program
	if ($mode == 2) {
		$tm = 4*$gc + 2*($len-$gc);
		$tm = sprintf("%.2f", $tm);
	}
	if ($mode == 3) {
		$oligotm = `which oligotm`;
		$tm = `oligotm -dv 2.5 -n 0.2 $seq`;
		$tm = sprintf("%.2f", $tm);
	}
	return $tm;
}
