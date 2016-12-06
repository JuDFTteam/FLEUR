package judft_tests;
@EXPORT="";

use jt;

#This runs a single test
sub execute_test($$$$){
    my $testdir=shift;
    my $exec=shift;
    my $mpi=shift;
    $workdir=shift;

    #is this test ok for the current configuration?
    my ($stages,$test_name)=
               test_applicable($testdir,$exec);
    if ($stages>0) {
	$executable="$mpi $exec";
	my $starttime=time();
	jt::initlog($workdir,$test_name,$executable);
    #prepare the workdir
	if (-r $workdir){
	    system("rm -f $workdir/*");
	}else{
	    die "Invalid workdir:$workdir" if (system("mkdir $workdir"));
	}
    
	my $old_dir=`pwd`;
	chomp $old_dir;
	chdir("tests/$testdir");
        #run all stages of the test
	for(my $stage=1;$stage<=$stages;$stage++){
	    print "Stage: $stage / $stages\n";
		do "test.run$stage";
	}   
	chdir($old_dir);
	my $time=time()-$starttime;
	jt::stoplog($time);
    }
    jt::testresult($workdir);
}    
    

sub test_applicable($$){
    my $testdir=shift;
    my $exec=shift;   


    #read description of test
    do "tests/$testdir/test.desc";

    #check if executable name starts with code name
    #if (!($exec=~/^\Q$test_code\E/i)){
    #	return (0,"");
    #}
    #test requirements

#    if ($exec=~/_SOC/){
#	return (0,"") if ($test_requirements{"SOC"}==0);
#    }else{
#	return (0,"") if ($test_requirements{"SOC"}==1);
#    }
#    if ($exec=~/_INVS/){
#	return (0,"") if ($test_requirements{"complex"}==1);
#    }
    if (!($exec=~/_MPI/)){
	return (0,"") if ($test_requirements{"MPI"}==1);
    }
    return ($test_stages,$test_name);
}
	
