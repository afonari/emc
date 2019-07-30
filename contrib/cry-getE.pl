#!/usr/bin/perl

#  VERSION 1.5
#
#
use strict;
use Data::Dumper;
use Getopt::Long;
    $Getopt::Long::ignorecase = 1;
use File::Copy;

use constant H2EV => 27.21138505; #
use constant ISS => 10000; # shrinking factor
# use constant NKPOINTS => 61; #

my ($opt_help, $opt_f9, $opt_nband);

GetOptions('f=s' => \$opt_f9, 'b=i' => \$opt_nband);

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
    
        $0 -f ../input.f9 -b 236

OUT
    exit;
}

# check for runprop09 in path
system("which runprop09 >& /dev/null") == 0 or die "runprop09 is not in \$PATH (`which` failed): $!";

my $band = $opt_nband;
print "Will use band = $band\n";

open( my $kpoints_fh, "<", "KPOINTS" ) || die "Can't open KPOINTS file: $!";

<$kpoints_fh>; # title
my $NKPOINTS = <$kpoints_fh> + 0; # make it int
<$kpoints_fh>; # Reciprocal

open( my $eigenval_fh, ">", "EIGENVAL" ) || die "$!\n";
print $eigenval_fh "$band\n";
print $eigenval_fh "LINE\n";
print $eigenval_fh "LINE\n";
print $eigenval_fh "LINE\n";
print $eigenval_fh "LINE\n";
print $eigenval_fh "1 1 1\n";


for(my $i = 1; $i <= $NKPOINTS; $i++)
{
    open( my $band_fh, ">", "input.d3" ) || die "Can't open input file: $!";

    print $band_fh "RDFMWF\n";
    print $band_fh "BAND\n";
    print $band_fh "For k-point (CART): TODO\n";
    print $band_fh "1 ".ISS." 2 $band $band 1 0\n";

    my @v1;
    ($v1[1], $v1[2], $v1[3]) = (<$kpoints_fh> =~ m/^\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);

    my @v2;
    if($i == $NKPOINTS)
    {   
        ($v2[1], $v2[2], $v2[3]) = (0.0, 0.0, 0.0)
    }
    else
    {
        ($v2[1], $v2[2], $v2[3]) = (<$kpoints_fh> =~ m/^\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
        $i++;
    }

    print $band_fh sprintf("%5d %5d %5d   %5d %5d %5d\n", $v1[1]*ISS, $v1[2]*ISS, $v1[3]*ISS, $v2[1]*ISS, $v2[2]*ISS, $v2[3]*ISS);
    print $band_fh "END\n";
    close($band_fh);

    print "\n\n\nrunning $i of ".$NKPOINTS.":\n\n";
    print `cat input.d3`;
    print "\n\n";
    copy($opt_f9,"input.f98") or die "Copy of $opt_f9 to ./input.f9 failed: $!";
    `runprop09 input input`;
    #`cat input.d3 > BAND-$i`;
    #`cat input_input_dat.BAND >> BAND-$i`;

    open( my $out_fh, "<", "input_input.outp" ) || die "Can't open file: $!";
    while(<$out_fh>)
    {
        if(/CARTESIAN COORD/){ print; }
        if(/POINTS/){ print; }
    }
    close($out_fh);
    
    open( my $bandout_fh, "<", "input_input_dat.BAND" ) || die "Can't open file: $!";
    my $c = 0;
    while(my $line = <$bandout_fh>)
    {
        if($line =~ /^\s+([-+\d\.E]+)\s+([-+\d\.E]+)\s+$/)
        {
            my $line1='';
            if( $c == 0){
                $line1 = sprintf("%5f %5f %5f\n", $v1[1], $v1[2], $v1[3]);
                $c++;
            }else{
                $line1 = sprintf("%5f %5f %5f\n", $v2[1], $v2[2], $v2[3])
            }
            print;
            # $line = <$bandout_fh>;
            print $eigenval_fh "\n";
            print $eigenval_fh $line1;
            print $eigenval_fh "1 ".sprintf("%15.10f\n", $2*H2EV); # to eVs
        }
    }
    close($bandout_fh);
}

print "\n\ndone.\n";
close($eigenval_fh);
