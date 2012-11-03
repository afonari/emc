#!/usr/bin/perl

use strict;
use Data::Dumper;
use Math::BigFloat; # for comparison of floats with the same precision set.
use Getopt::Long;
	$Getopt::Long::ignorecase = 1;

use constant PI => 3.14159265358979324;
use constant A2B => 1.88972613289; # ANGSTROM_TO_BOHR
use constant B2A => 0.529177249; #

my (@opt_kp, $opt_help, $opt_input);
our ($dk, @et, %ek);

GetOptions('kp=f{3}' => \@opt_kp, 'help' => \$opt_help);

if(scalar(@opt_kp)!=3)
{
    print "Please input k point (see below) !!\n";
    $opt_help = 1;
}

if($opt_help)
{
    print <<OUT;
    
    Run as:
    
        $0 -k 0.5 0.0 0.5

OUT
    exit;
}

open( my $inp_fh, "<", "inp" ) || die "Can't open inp file: $!";
<$inp_fh>;
<$inp_fh>;
<$inp_fh>;
my $prg = <$inp_fh>;
$prg = trim($prg);

open( my $outcar_fh, "<", "OUTCAR" ) || die "Can't open OUTCAR: $!";

my (@f, @g);
while(<$outcar_fh>)
{
    if($prg eq 'C')
    {
        if(/DIRECT LATTICE VECTORS COMPON. \(A.U.\)\s+RECIP. LATTICE VECTORS COMPON. \(A.U.\)/)
        {
            print "Will treat OUTCAR as CRYSTAL output.\n";
            $_ = <$outcar_fh>; # next line
            ($f[1], $f[2], $f[3], $g[1], $g[2], $g[3]) = (<$outcar_fh> =~ m/\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
            ($f[4], $f[5], $f[6], $g[4], $g[5], $g[6]) = (<$outcar_fh> =~ m/\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
            ($f[7], $f[8], $f[9], $g[7], $g[8], $g[9]) = (<$outcar_fh> =~ m/\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
            last;
        }
    }
    elsif($prg eq 'V')
    {
        if(/direct lattice vectors/)
        {
            print "Will treat OUTCAR as VASP output.\n";
            ($f[1], $f[2], $f[3], $g[1], $g[2], $g[3]) = (<$outcar_fh> =~ m/\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
            ($f[4], $f[5], $f[6], $g[4], $g[5], $g[6]) = (<$outcar_fh> =~ m/\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
            ($f[7], $f[8], $f[9], $g[7], $g[8], $g[9]) = (<$outcar_fh> =~ m/\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
            last;
        }
    }    
}
close $outcar_fh;

print <<OUT;

f:
  $f[1], $f[2], $f[3],
  $f[4], $f[5], $f[6],
  $f[7], $f[8], $f[9]

g:
  $g[1], $g[2], $g[3],
  $g[4], $g[5], $g[6],
  $g[7], $g[8], $g[9]

OUT

@f = map { $_ * 1.0/(2*PI) } @f;

my @kp = (undef, $opt_kp[0], $opt_kp[1], $opt_kp[2]);
print "from input: ".sprintf("%10.7f %10.7f %10.7f\n", $kp[1], $kp[2], $kp[3]);

my @kpr = m3by1( t3by3(@g), @kp );
print "T(g)*kp: ".sprintf("%10.7f %10.7f %10.7f\n", $kpr[1], $kpr[2], $kpr[3]);


@kpr = m3by1( @g, @kp ); # to reciprocal cartesian space
@kpr = n1(@kpr);
print "Norm(g*kp): ".sprintf("%10.7f %10.7f %10.7f\n", $kpr[1], $kpr[2], $kpr[3]);


# sum two three-component vectors
sub s1w1(@)
{
    my (@m, @n);
    
    @m[0..3] = @_[0..3];
    @n[0..3] = @_[4..7];
    
    return (undef, $m[1]+$n[1], $m[2]+$n[2], $m[3]+$n[3]);
}

sub n1(@)
{
    my (@m);
    
    @m[0..3] = @_[0..3];
    my $n = sqrt($m[1]**2 + $m[2]**2 + $m[3]**2);
    return (undef, $m[1]/$n, $m[2]/$n, $m[3]/$n);
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

sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}