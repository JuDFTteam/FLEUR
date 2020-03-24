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
   te.log_info("GaAsHybridPBE0 test")
   copyfile(f"{test_loc}/files/inp.xml", f"{te.workdir}/inp.xml")

   # special for this hybrid test:
   te.nprocs = 3

   if(te.parallel):
      te.run(["mpirun", "-n", f"{te.nprocs}", "--allow-run-as-root", te.binary])
   else:
      te.run([te.binary])


   te.check_value_outfile("HF total energy=", "htr", [-4205.2193168, -4204.9174214647], 0.000001)

   # only check the last bandgap
   exp_bandgap = 38 * [None]
   exp_bandgap[37] = 0.09105
   te.check_value_outfile("bandgap                     :", "htr", exp_bandgap, 0.0001)

except Exception as e:
   te.log_error(f"Error: {e}")
   sys.exit(1)