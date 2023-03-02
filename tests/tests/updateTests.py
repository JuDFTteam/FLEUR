import re 
import shutil
with open("Testing/pytest_session.stdout","r") as f:
    reg_file=re.compile("out.xml file:(\S+) against values in \.(\S+)")
    reg_testname=re.compile("________ (\S+) _________")
    reg_failed=re.compile("Check failed for ")
    out_test=""
    out_comp=""
    for s in f:
        m=reg_testname.search(s)
        if m:
            testdir=m.groups()[0]
        m=reg_file.search(s)
        if m:
            out_test=m.groups()[0]
            out_comp="../tests"+m.groups()[1]            
            m=re.search("(.*)work(.*)",out_test)
            out_test=m.groups()[0]+"failed_test_results/"+testdir+m.groups()[1]
        if  reg_failed.search(s):   
            print("Copy:",out_test,out_comp)
            shutil.copy(out_test,out_comp)
        
