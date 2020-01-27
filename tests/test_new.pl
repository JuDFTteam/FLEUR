#!/usr/bin/perl
# Driver for cmake testing of FLEUR
# See README.tests for details

use strict;
my $dir;
BEGIN{
$dir=$0;
$dir=~s/test_new.pl//;
}
use lib "$dir/scripts";

use judft_tests;

my $testdir=shift;
my $executable=shift;

#check MPI environment
my $mpi=shift;
if ($mpi){
    if ($ENV{"juDFT_MPI"}) {
	    $mpi=$ENV{"juDFT_MPI"};	
    }  
}
if ($executable=~/_MPI/){
   if (!$mpi){
     #Try default mpi setting if none was given
      $mpi="mpirun -np 2 "
   }
}else{
    #no mpi executable...
    $mpi='';
}
print "MPI:$mpi\n";

my $workdir="$ENV{PWD}/Testing/work";
system("rm $workdir/*");
chdir($dir);

judft_tests::execute_test($testdir,$executable,$mpi,$workdir);
