import sys
import argparse
import os
import logging
import shutil
import subprocess

class TestEnv:
   binary = ""
   testdir = ""
   workdir = ""
   command = ""
   parallel = False
   nprocs = -1
   errors = 0

   def __init__(self):
      parser = argparse.ArgumentParser(description='get test-dir and bin-dir')
      parser.add_argument("--bindir", type=str, nargs=1, default='', required=False, help="location of fleur executable")
      parser.add_argument("--command", type=str, nargs=1, default='', required=False, help="command, which to execute instead of fleur")
      parser.add_argument("--testdir", type=str, nargs=1, required=True, help="where to execute tests")
      parser.add_argument("--nprocs", type=int, nargs='?', default=2, help="number of parallel mpi-processes")
      args = parser.parse_args()

      self.setup_logger(args)
      self.setup_env(args)
      if args.bindir:
         self.find_binary(args)
      elif args.command:
         self.command = args.command[0]
      else:
         logging.warning("No executable/command specified")
         sys.exit(1)
      self.nprocs = args.nprocs

   def find_binary(self, args):      
      fleur_dir = args.bindir[0]
      if(fleur_dir[-1] == "/"):
         fleur_dir = fleur_dir[:-1]

      if(os.path.isfile(f"{fleur_dir}/fleur_MPI")):
         self.binary = f"{fleur_dir}/fleur_MPI"
         logging.info(f"Use {self.binary} as executable")
         self.parallel = True

      elif(os.path.isfile(f"{fleur_dir}/fleur")):
         self.binary = f"{fleur_dir}/fleur"
         logging.info(f"Use {self.binary} as executable") 
         self.parallel = False

      elif(os.path.isfile(f"{fleur_dir}/inpgen2")):
         self.binary = f"{fleur_dir}/inpgen2"
         logging.info(f"Use {self.binary} as executable")
         self.parallel = False

      else:
         logging.warning("Can not find any executables")
         sys.exit(1)
   
   def setup_logger(self,args):
      if not os.path.isdir(args.testdir[0]):
         os.makedirs(args.testdir[0])
      logging.basicConfig(filename=f"{args.testdir[0]}/test.log",level=logging.DEBUG, format='%(asctime)s %(message)s')
      logging.info("###############################################################")

   def setup_env(self, args):
      self.testdir = args.testdir[0]
      if(self.testdir[-1] == "/"):
         self.testdir = self.testdir[:-1]

      self.workdir = f"{self.testdir}/work"
      if(os.path.isdir(self.workdir)):
         shutil.rmtree(self.workdir)
      os.makedirs(self.workdir)
      os.chdir(self.workdir)
   
   def log_info(self, text):
      logging.info(text)
   
   def log_warn(self, text):
      logging.warn(text)
   
   def log_error(self, text):
      logging.error(text)

   def run(self,arg_list):
      if not self.binary:
         self.log_error('No fleur executable specified')
      OMP_NUM_THREADS = os.environ['OMP_NUM_THREADS'] if 'OMP_NUM_THREADS' in os.environ else ''
      self.log_info(f"Start running command:\nOMP_NUM_THREADS={OMP_NUM_THREADS}\n{arg_list}")
      with open(f"{self.workdir}/stdout", "w") as f_stdout:
         with open(f"{self.workdir}/stderr", "w") as f_stderr:
            subprocess.run(arg_list + ["-no_send"], stdout=f_stdout, stderr=f_stderr, check=True)
      self.log_info("Finished running")
      try:
         with open(f"{self.workdir}/out.xml", "r") as outxml:
            for line in outxml.readlines():
               if "ERROR Message" in line:
                  self.log_info(line.strip())
      except OSError:
         pass

   
   def check_value_outfile(self, before_str, after_str, expected, delta):
      exp_idx = 0
      with open(f"{self.workdir}/out", "r") as f:
         found = False
         for line in f.readlines():
            if(before_str in line):
               value_string = line.split(before_str)[-1]
               value_string = value_string.split(after_str)
               # remove empty strings
               value_string = [i for i in value_string if i is not ""][0]
               value = float(value_string)
               
               if(expected[exp_idx] is not None):
                  if(abs(value - expected[exp_idx]) < delta):
                     self.log_info(f"PASSED [{exp_idx}]: {before_str} found: {value} expected: {expected[exp_idx]}")
                  else:
                     self.log_info(f"FAILED [{exp_idx}]: {before_str} found: {value} expected: {expected[exp_idx]}")
                     self.errors += 1
               exp_idx += 1
               found = True
      if(len(expected) != exp_idx):
         self.log_error("Number of expected values disagree with found values.")
         self.log_error(f"before_str: '{before_str}''")
         sys.exit(1)
      if(not found):
         self.log_error(f"{before_str} not found in file")
         sys.exit(1)

      if(self.errors != 0):
         sys.exit(1)

