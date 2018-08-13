# to be executed in the fleur/ directory

from glob import glob

files = []
files.extend(glob("**/*.f90"))
files.extend(glob("**/*.F90"))

cmake_files = glob("**/CMakeLists.txt")

def rm_comments(line):
  comm_pos = line.find("#")
  if comm_pos == -1:
    return line
  else:
    return line[:comm_pos]


def in_cmake(cmake_files, file):
  for cmf in cmake_files:
    with open(cmf) as f:
      lines = f.readlines()
      for l in lines:
        if file in rm_comments(l):
          return True
  return False



for f in files:
  if(not(in_cmake(cmake_files, f))):
    print(f)