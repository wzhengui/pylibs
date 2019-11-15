#!/usr/bin/perl
#
# Convert SMS 2dm grid to .gr3 format, simplified. Assumes no extra nodes.
#

if ( @ARGV != 2 ) {
    print "Usage: $0 [SMS grid file] [.gr3 file]\n";
    exit(1);
}

$file    = $ARGV[0];
$outfile = $ARGV[1];

open( IN,  "<$file" )    || die "Can't open $file\n";
open( OUT, ">$outfile" ) || die "Can't open $outfile\n";

my @lines = <IN>;
close(IN);

my $ne = 1;
my $np = 1;
my @nes = ();
my @nps = ();

for (my $i=0;$i<@lines; $i++) {
    if ($lines[$i] =~ /^E3T/) {
       my ($junk, $junk, $n1, $n2, $n3, $junk) = split(' ', $lines[$i]);
       $nes[$ne - 1] = "$ne 3 $n1 $n2 $n3\n";
       $ne++;
    } elsif ($lines[$i] =~ /^E4Q/){
       my ($junk, $junk, $n1, $n2, $n3,$n4, $junk) = split(' ', $lines[$i]);
       $nes[$ne - 1] = "$ne 4 $n1 $n2 $n3 $n4\n";
       $ne++;
    } elsif ($lines[$i] =~ /^ND/) {
       my ($junk, $num, $x, $y, $d) = split(' ', $lines[$i]);
       $nps[$np - 1] = "$np $x $y $d\n";
       $np++;
    }
}

$ne--;
$np--;

print OUT "2dm2gr3\n$ne $np\n";;
print OUT "@nps";
print OUT "@nes";
close(OUT);
