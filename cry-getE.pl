#!/usr/bin/perl

#  VERSION 1.0 (or 1.1.5 of the old effective mass script)
#
#
use strict;
use Data::Dumper;
use Math::BigFloat; # for comparison of floats with the same precision set.
use Getopt::Long;
	$Getopt::Long::ignorecase = 1;

use constant PI => 3.14159265358979324;
use constant A2B => 1.88972613289; # ANGSTROM_TO_BOHR
use constant B2A => 0.529177249; #
use constant ISS => 1000; # shrinking factor
use constant NKPOINTS => 61; #

my ($opt_help, $opt_input, $opt_f9, $opt_nband);

GetOptions('i=s' => \$opt_input, 'f=s' => \$opt_f9, 'b=i' => \$opt_nband);

if(!$opt_input)
{
    print "Please input SCF.out\n";
    $opt_help = 1;
}

if(!$opt_f9)
{
    print "Please input SCF.f9\n";
    $opt_help = 1;
}

if(!$opt_nband)
{
    print "Please input band number\n";
    $opt_help = 1;
}

if($opt_help)
{
    print <<OUT;
    
    Run in a child directory of a (presumably) SCF run.
    Example:
    
        $0 -i ../input.out -f ../input.f9 -b 236

OUT
    exit;
}

# check for runprop09 in path
system("which runprop09 >& /dev/null") == 0 or die "runprop09 is not in \$PATH (which call failed): $?";

open( my $outcar_fh, "<", $opt_input ) || die "Can't open $opt_input: $!";

my (@f, @g, $vb, $cb);
while(<$outcar_fh>)
{
    if(/DIRECT LATTICE VECTORS COMPON. \(A.U.\)\s+RECIP. LATTICE VECTORS COMPON. \(A.U.\)/)
    {
        $_ = <$outcar_fh>; # next line

        # hate CRYSTAL for that! ;)
        ($f[1], $f[2], $f[3], $g[1], $g[2], $g[3]) = (<$outcar_fh> =~ m/\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
        ($f[4], $f[5], $f[6], $g[4], $g[5], $g[6]) = (<$outcar_fh> =~ m/\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
        ($f[7], $f[8], $f[9], $g[7], $g[8], $g[9]) = (<$outcar_fh> =~ m/\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
    }
    if(/N. OF ELECTRONS PER CELL\s*(\d+)/)
    {
        $vb = $1/2;
        $cb = $vb+1;

    }
}
close $outcar_fh;

my $band = $opt_nband;
print "Will use band = $band\n";


# @f = map { $_ * 1.0/(2*PI) } @f;

print <<OUT;

direct lattice vectors matrix (from $opt_input):
  $f[1], $f[2], $f[3],
  $f[4], $f[5], $f[6],
  $f[7], $f[8], $f[9]

reciprocal lattice vectors matrix:
  $g[1], $g[2], $g[3],
  $g[4], $g[5], $g[6],
  $g[7], $g[8], $g[9]

OUT

open( my $kpoints_fh, "<", "KPOINTS" ) || die "Can't open KPOINTS file: $!";

<$kpoints_fh>; # title
<$kpoints_fh>; # NKPOINTS == 61
<$kpoints_fh>; # Cartesian

open( my $eigenval_fh, ">", "EIGENVAL" ) || die "$!\n";
print $eigenval_fh "$band\n";
print $eigenval_fh "LINE\n";
print $eigenval_fh "LINE\n";
print $eigenval_fh "LINE\n";
print $eigenval_fh "LINE\n";
print $eigenval_fh "1 1 1\n";


for(my $i = 1; $i <= NKPOINTS; $i++)
{
    open( my $band_fh, ">", "input.d3" ) || die "Can't open input file: $!";

    print $band_fh "BAND\n";
    print $band_fh "For k-point (CART): TODO\n";
    print $band_fh "1 ".ISS." 2 $band $band 1 0\n";

    my @v1;
    ($v1[1], $v1[2], $v1[3]) = (<$kpoints_fh> =~ m/^\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);

    my @v2;
    if($i == NKPOINTS)
    {   
        ($v2[1], $v2[2], $v2[3]) = (0.0, 0.0, 0.0)
    }
    else
    {
        ($v2[1], $v2[2], $v2[3]) = (<$kpoints_fh> =~ m/^\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
        $i++;
    }

    my @kp1 = m3by1( t3by3(@f), @v1 );
    my @kp2 = m3by1( t3by3(@f), @v2 );
    print $band_fh sprintf("%5d %5d %5d   %5d %5d %5d\n", $kp1[1]*ISS, $kp1[2]*ISS, $kp1[3]*ISS, $kp2[1]*ISS, $kp2[2]*ISS, $kp2[3]*ISS);
    print $band_fh "END\n";
    close($band_fh);

    print "\n\n\nrunning $i of ".NKPOINTS.":\n\n";
    print `cat input.d3`;
    print "\n\n";
    `cp $opt_f9 ./`;
    `runprop09 input input`;
    #`cat input.d3 > BAND-$i`;
    #`cat input_input_dat.BAND >> BAND-$i`;

    open( my $out_fh, "<", "input_input.outp" ) || die "Can't open file: $!";
    while(<$out_fh>)
    {
        if(/CARTESIAN COORD/)
        {
            print;
        }

        if(/POINTS/)
        {
            print;
        }
    }
    close($out_fh);
    
    open( my $bandout_fh, "<", "input_input_dat.BAND" ) || die "Can't open file: $!";
    while(my $line = <$bandout_fh>)
    {
        if($line =~ /^\s+([-+\d\.E]+)\s+([-+\d\.E]+)\s+$/)
        {
            print;
            # $line = <$bandout_fh>;
            print $eigenval_fh "\n";
            print $eigenval_fh "LINE";
            print $eigenval_fh "\n";
            print $eigenval_fh "1 ".sprintf("%15.10f\n", $2*27.211399); # to eVs like VASP
        }
    }
    close($bandout_fh);
}

print "done.\n";
close($eigenval_fh);

# sum two three-component vectors
sub s1w1(@)
{
    my (@m, @n);
    
    @m[0..3] = @_[0..3];
    @n[0..3] = @_[4..7];
    
    return (undef, $m[1]+$n[1], $m[2]+$n[2], $m[3]+$n[3]);
}

sub m3by1(@)
{
    my (@m, @n, @r);
    
    @m[0..9] = @_[0..9];
    @n[0..3] = @_[10..13];
    
    $r[1] = $m[1]*$n[1] + $m[2]*$n[2] + $m[3]*$n[3];
    $r[2] = $m[4]*$n[1] + $m[5]*$n[2] + $m[6]*$n[3];
    $r[3] = $m[7]*$n[1] + $m[8]*$n[2] + $m[9]*$n[3];
    
    return @r;
}

sub t3by3(@)
{
    my (@m, @t);
    
    @m[0..9] = @_[0..9];
    $t[1] = $m[1]; $t[2] = $m[4]; $t[3] = $m[7];
    $t[4] = $m[2]; $t[5] = $m[5]; $t[6] = $m[8];
    $t[7] = $m[3]; $t[8] = $m[6]; $t[9] = $m[9];
       
    return @t;
}
