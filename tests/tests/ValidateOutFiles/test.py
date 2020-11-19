#!/usr/bin/env python3
# load libtest
import pathlib
test_loc = pathlib.Path(__file__).parent.absolute()
import sys
sys.path.insert(0, f"{test_loc}/../../")
from libtest import TestEnv
import os
import subprocess
from shutil import which

if which('xmllint') is None:
   te.log_error('Error: xmllint command is not available')
   sys.exit(1)

te = TestEnv()
te.log_info("OutputSchema Test")
with open(f"{te.workdir}/xmllint_out", "w") as f_stderr:
   test_dir = os.path.abspath(os.path.join(te.workdir,'./../..'))
   te.log_info(f"Test directory: {test_dir}")
   for root, dirs, files in os.walk(test_dir):
      for file in files:
         if file == 'out.xml':
            file_path = os.path.join(root,file)
            schema_path = os.path.join(root,'FleurOutputSchema.xsd')
            te.log_info(f"Testing {file_path}")
            if not os.path.isfile(schema_path):
               te.log_info("No OutputSchema found")
               continue
            with open(f"{te.workdir}/last_stdout", "w") as f_stdout:
               try:
                  subprocess.run(['xmllint', '--schema', f'{schema_path}', f'{file_path}'],
                                 stdout=f_stdout, stderr=f_stderr, check=True)
               except Exception as e:
                  test_name = root.replace(test_dir,'')
                  test_name = test_name.replace('work','')
                  test_name = test_name.replace('/','')
                  te.log_error(f"{test_name} failed to validate")
                  te.errors += 1

sys.exit(te.errors)

