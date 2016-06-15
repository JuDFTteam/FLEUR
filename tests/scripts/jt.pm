package jt;
@EXPORT="";

use POSIX;
use strict;
sub initlog($$$){
    my $workdir=shift;
    my $test_name=shift;
    my $config_name=shift;
    if (-r "$workdir/../test.log"){
	system("cat $workdir/../test.log >>$workdir/../test.oldlogs");
    }
    open(LOG,">$workdir/../test.log");
    print LOG "Configuration: $config_name\n";
    print LOG "*************************************\n";
    print LOG POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime);
    print LOG "Starting\n";
    print  "Configuration: $config_name\n";
    print  "Test: $test_name\n";
    print  "Workdir: $workdir\n";     
}

sub stoplog($){
    my $time=shift;
    print LOG "Finished after: $time\n";
    close(LOG);
}

sub copyfile($$){
    my $from=shift;
    my $to=shift;
    use POSIX;
    print LOG POSIX::strftime("%m/%d/%Y %H:%M:%S--", localtime);
    print LOG "Copying $from -> $to:";

    system("cp $from $to");
    my $res=system("diff -q $from $to");
    
    if ($res==0) {print LOG "Done\n";}
       else {print LOG "Failed\n";}
}
sub movefile($$){
    my $from=shift;
    my $to=shift;
    use POSIX;
    print LOG POSIX::strftime("%m/%d/%Y %H:%M:%S--", localtime);
    print LOG "Moving $from -> $to:";

    system("mv $from $to");
    print LOG "Done\n";
}
sub deletefile($){
    my $from=shift;
    use POSIX;
    print LOG POSIX::strftime("%m/%d/%Y %H:%M:%S--", localtime);
    print LOG "Deleting $from";

    system("rm $from ");
    print LOG "Done\n";
}

sub testrun($$){
    my $ex=shift;
    my $dir=shift;

    print LOG POSIX::strftime("%m/%d/%Y %H:%M:%S--", localtime);
    print LOG "Running $ex:";
    
    if (system("cd $dir;$ex")==0){
	print LOG "Done\n";}
       else {
	   print LOG "Failed\n";}

    print LOG POSIX::strftime("%m/%d/%Y %H:%M:%S--", localtime);
    print LOG "Finished execution\n";
}

sub test_fileexists($){
    my $file=shift;
    print LOG POSIX::strftime("%m/%d/%Y %H:%M:%S--", localtime);
    print LOG "Testing for $file:";
    
    if (-r $file){
	print LOG "Exists\n";
	return 0;
    }else{
	print LOG "Not found\n";
	return 1;
    }
}

sub test_grepexists($$){
    my $file=shift;
    my $text=shift;
    print LOG POSIX::strftime("%m/%d/%Y %H:%M:%S--", localtime);
    print LOG "Grep for $text in $file:";

    if (system("grep -q \"$text\" $file")==0){
    print LOG "Found\n";
	return 0;
    }else{
	print LOG "Not found\n";
	return 1;
    }
}

sub test_grepnumber($$$$$){
    my $file=shift;
    my $grepfor=shift;
    my $reg=shift;
    my $value=shift;
    my $tol=shift;

    print LOG POSIX::strftime("%m/%d/%Y %H:%M:%S--", localtime);
    print LOG "Grep for $grepfor in $file:";

    my $l=`grep \"$grepfor\" $file`;

    
    $l=~m/$reg/s;

    
    print LOG "$1 == $value:";

    if (abs($1-$value)<$tol){
	print LOG "ok\n";
	return 0;
    }else{
	print LOG "failed\n";
	return 1;
    }
}

       


sub stageresult($$$){
    my $workdir=shift;
    my $result=shift;
    my $stage=shift;
    
    system("rm -f $workdir/test_$stage.*");
    if ($result==0){
	system("touch $workdir/test_$stage.ok");
	print "Stage $stage passed\n";
    }else{
	system("touch $workdir/test_$stage.failed");
	print "Stage $stage failed\n";
    }
    
}

sub testresult($){
    my $workdir=shift;
    my $result;

    system("rm -f $workdir/test.ok $workdir/test.failed");
    $result=1 if (system("ls $workdir/test_?.failed 2>/dev/null")==0);

    if ($result==0){
	system("touch $workdir/test.ok");
	print "Test passed\n";
	exit 0;

    }else{
	system("touch $workdir/test.failed");
	print  "Test failed\n";
	die "Test failed";
    }
}

    
