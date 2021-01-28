#!/usr/bin/perl
#Convert gredit grid file to SMS 2dm file
if (@ARGV !=2) {
	die "usage: $0 infile outfile\n";
}
$file = $ARGV[0];

$outfile = $ARGV[1];
open(IN,$file);
@lines = <IN>;
close(IN);
open(OUT,">$outfile");
print OUT "MESH2D\n";
chomp $lines[1];
$lines[1]=~s/^\s+//;
$lines[1]=~s/\s+/ /g;
($e,$n)=split(" ",$lines[1]);
$starte = $n+2;
print "$lines[1]\n$e $n $starte\n";
for ($i = $starte; $i<$n+$e+2; $i++){
chomp $lines[$i];
$lines[$i]=~s/^\s+//;
$lines[$i]=~s/\s+/ /g;
	($elemn,$elem34,$e1,$e2,$e3,$e4)=split(" ",$lines[$i]);
	if ($elem34 == 4) {print OUT "E$elem34"."Q $elemn $e1 $e2 $e3 $e4 1\n";}
	elsif ($elem34 == 3) {print OUT "E$elem34"."T $elemn $e1 $e2 $e3 1\n";}
}
for ($i = 2; $i<$starte; $i++){
chomp $lines[$i];
$lines[$i]=~s/^\s+//;
$lines[$i]=~s/\s+/ /g;
	print OUT "ND $lines[$i]\n";
}
