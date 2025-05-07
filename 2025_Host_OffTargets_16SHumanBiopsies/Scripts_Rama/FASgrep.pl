#!/usr/bin/perl
if (!$ARGV[0]) {
	print "USAGE: FASgrep.pl [-i] [-v] [-c] [-s] [-f ][file] key\n\n";
	print "       -i   : case insensitive\n";
	print "       -v   : negate match\n";
	print "       -c   : just count, not print\n";
	print "       -s   : search in sequence\n";
	print "       -m   : print only the match, not the whole sequence\n";
	print "       -f   : use keys in file for exact matching\n";
	print "       -p   : input sequences are proteins, not DNA\n";
	print "       file : a file to scan into\n";
	print "       key  : the search key\n";
	print "       -d   : debug mode\n\n";
	exit;
}

$inseq = 0;
$case = 1;
$nega = 0;
$cont = 0;
$isprot = 0;

$count = 0;
foreach $o (@ARGV) {
	#the "-" marks options
	if ($o =~ /^-/) {
		$case = 0 if ($o =~ /-i/);
		$nega = 1 if ($o =~ /-v/);
		$cont = 1 if ($o =~ /-c/);
		$inseq = 1 if ($o =~ /-s/);
		$isprot = 1 if ($o =~ /-p/);
		$match = 1 if ($o =~ /-m/);
		$debug = 1 if ($o =~ /-d/);
		&degen if ($o =~ /-s$/);
		if ($o =~ /-f$/) {
			shift @ARGV;
			#this is a file containing search patterns
			$key = shift @ARGV;
			next;
		}
	}
	#without "-"
	else {
		#this is the sequence file
		open(STDIN,$o) if (-e $o && $o !~ /\.gz$/);
		print STDERR "gzip input detected\n" if (-e $o && $o =~ /\.gz$/);
 		open(STDIN,"zcat $o|") if (-e $o && $o =~ /\.gz$/);
		#this is the search pattern(s)
		push(@patt,$o) if (!-e $o);
	}
}
#print STDERR "sequences are prot\n" if ($isprot);
if ($key) {
	print STDERR "Loading keys from $key...";
	open(IN,$key);
	while($line = <IN>) {
		$line =~ s/[\n\r]//g;
		next if ($line !~ /\w/);
		$ok{$line} = 1;
	}
	close IN;
	$arg = join "|",keys %ok;
	print STDERR "done.\n";
}

if (scalar @patt == 1) {
	$arg = $patt[0];
}

if ($inseq) {
#	print STDERR "degenerating...\n";
	foreach $p (0..$#patt) {
#		print $patt[$p],"--";
		$patt[$p] =~ s/(\w)/$deg{$1}/gi if (!$isprot);
#		print $patt[$p],"\n";
	}
	$arg = join ".+?",@patt;
}

$arg = lc($arg) if ($case == 0);
print $arg,"\n" if ($debug);
#exit;

while($l = <STDIN>) {
	chomp $l;
	if ($l =~ />/) {
		if ($seq) {#so at the second ">" 
			scan($seq);
		}
		$name = $l;
		$seq = '';
	} else {	
		$seq .= $l;
#		print STDERR ">$name\n$seq\n";
#		sleep 1;
	}
}
scan($seq);

$print = 0;
$q = $name if ($inseq == 0);
$q = $seq if ($inseq == 1);
$q = lc($q) if ($case == 0);
$print = 1-$nega if ($q =~ /$arg/);
$print = 0-$nega if ($q !~ /$arg/);
$count++ if ($print != 0);
$seq = $1 if ($print && $match && !$nega);
print $name,"\n",$seq,"\n" if ($print != 0 && $cont == 0);

print $count,"\n" if ($cont == 1);

sub scan {
	my $seq = shift;
	$print = 0;
	$q = $name if ($inseq == 0);
	$q = $seq if ($inseq == 1);
	$q = lc($q) if ($case == 0);
	#$print = match($q,$arg);
	if ($q =~ /($arg)/) {
		$print = 1-$nega;
		$seq = $1 if ($print && $match && !$nega);
	} else {
		$print = 0-$nega;
	}
	$count++ if ($print != 0);
#	$print = 1 if ($ok{$name});
	#print STDERR "Q: $arg | $q : $print $count\n";
#	sleep 10;
	print $name,"\n",$seq,"\n" if ($print != 0 && $cont == 0);
#	<STDIN>;
}

#sub match {
#	my $q = shift;
#	my $arg = shift;
#	$print = 1-$nega if ($q =~ /$arg/);
#	$print = 0-$nega if ($q !~ /$arg/);
#	print STDERR "$q - $arg -> $print";
#	<STDIN>;
#	return $print;
#}

sub degen {
(%deg) = (
"A" => "A",
"T" => "T",
"C" => "C",
"G" => "G",
"R" => "[AG]",
"Y" => "[CT]",
"S" => "[CG]",
"W" => "[AT]",
"K" => "[GT]",
"M" => "[AC]",
"B" => "[CGT]",
"D" => "[AGT]",
"H" => "[ACT]",
"V" => "[ACG]",
"N" => "[ACGT]",

"a" => "a",
"t" => "t",
"c" => "c",
"g" => "g",
"r" => "[ag]",
"y" => "[ct]",
"s" => "[cg]",
"w" => "[at]",
"k" => "[gt]",
"m" => "[ac]",
"b" => "[cgt]",
"d" => "[agt]",
"h" => "[act]",
"v" => "[acg]",
"n" => "[acgt]"
);
}
