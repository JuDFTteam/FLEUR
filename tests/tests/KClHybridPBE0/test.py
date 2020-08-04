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
   te.log_info("KClHybridPBE0 test")
   copyfile(f"{test_loc}/files/inp.xml", f"{te.workdir}/inp.xml")

   # special for this hybrid test:
   te.nprocs = 3


   if(te.parallel):
      te.run(["mpiexec", "-n", f"{te.nprocs}", te.binary, "-trace"])
   else:
      te.run([te.binary, "-trace"])


   te.check_value_outfile("HF total energy=", "htr", [-1063.8587731477, -1063.8383730939], 0.000001)

   # only check the last bandgap
   exp_bandgap = 28 * [None]
   exp_bandgap[27] = 0.2771
   te.check_value_outfile("bandgap                     :", "htr", exp_bandgap, 0.0001)

   sys.exit(te.errors)

except Exception as e:
   te.log_error(f"Error: {e}")
   sys.exit(1)
