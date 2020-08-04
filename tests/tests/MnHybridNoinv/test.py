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
   te.log_info("MnHybridNoinv test")
   copyfile(f"{test_loc}/files/inp.xml", f"{te.workdir}/inp.xml")

   # special for this hybrid test:
   te.nprocs = 3

   if(te.parallel):
      te.run(["mpiexec", "-n", f"{te.nprocs}", , te.binary])
   else:
      te.run([te.binary])


   te.check_value_outfile("HF total energy=", "htr", [-2317.16697943, -2317.1756358345], 0.000001)

   exp_mm1 = 32*[None]
   exp_mm1[-1] = 2.47449
   te.check_value_outfile("--> mm       1", " ", exp_mm1, 0.00001)

   exp_mm2 = 32*[None]
   exp_mm2[-1] =-2.47449
   te.check_value_outfile("--> mm       2", " ", exp_mm2, 0.00001)

   # only check the last bandgap
   exp_bandgap = 32 * [None]
   exp_bandgap[-1] = 0.087531
   te.check_value_outfile("bandgap                     :", "htr", exp_bandgap, 0.0001)

   sys.exit(te.errors)

except Exception as e:
   te.log_error(f"Error: {e}")
   sys.exit(1)
