#!/usr/bin/env python3
# load libtest
import pathlib
test_loc = pathlib.Path(__file__).parent.absolute()
import sys
sys.path.insert(0, f"{test_loc}/../../")
from libtest import TestEnv
import os
import subprocess

te = TestEnv()
te.log_info("OutputSchema Test")
with open(f"{te.workdir}/xmllintErrors", "w") as f_stderr:
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
               te.errors += 1
               continue
            with open(f"{te.workdir}/last_xmllintOut", "w") as f_stdout:
               arg_list = [te.command, '--schema', f'{schema_path}', f'{file_path}']
               te.log_info(f"Running command: {' '.join(arg_list)}")
               try:
                  subprocess.run(arg_list, stdout=f_stdout, stderr=f_stderr, check=True)
               except Exception as e:
                  te.log_error(f'Error: {e}')
                  test_name = root.replace(test_dir,'')
                  test_name = test_name.replace('work','')
                  test_name = test_name.replace('/','')
                  te.log_error(f"{test_name} failed to validate")
                  te.errors += 1

sys.exit(te.errors)

