#!/usr/bin/perl
# Driver for cmake testing of FLEUR
# See README.tests for details

use strict;
my $dir;
BEGIN{
$dir=$0;
$dir=~s/test.pl//;
}
use lib "$dir/scripts";

use judft_tests;

my $testdir=shift;
my $executable=shift;
my $mpi=shift;
if (!$mpi){
    $mpi=$ENV{"juDFT_MPI"};
}
print "MPI:$mpi\n";
my $workdir="$ENV{PWD}/Testing/work"; 
chdir($dir);

judft_tests::execute_test($testdir,$executable,$mpi,$workdir);
