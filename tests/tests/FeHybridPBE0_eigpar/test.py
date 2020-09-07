#!/usr/bin/env python3
# load libtest
import pathlib
test_loc = pathlib.Path(__file__).parent.absolute()
import sys
sys.path.insert(0, f"{test_loc}/../../")
from libtest import TestEnv
from shutil import copyfile

try:
   te = TestEnv()
   te.log_info("FeHybridPBE0 test")
   copyfile(f"{test_loc}/files/inp.xml", f"{te.workdir}/inp.xml")

   # special for this hybrid test:
   te.nprocs = 6

   if(te.parallel):
      te.run(["mpirun", "-n", f"{te.nprocs}", "--allow-run-as-root", "--mca", "btl", "vader,self" , te.binary, "-trace"])
   else:
      te.run([te.binary, "-trace"])


   te.check_value_outfile("HF total energy=", "htr", [-1272.72894, -1272.7324763242], 0.000001)

   exp_mm1 = 27*[None]
   exp_mm1[14] = 3.40037
   exp_mm1[-1] = 3.42075
   te.check_value_outfile("--> mm       1", " ", exp_mm1, 0.00001)

   sys.exit(te.errors)

except Exception as e:
   te.log_error(f"Error: {e}")
   sys.exit(1)
