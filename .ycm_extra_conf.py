import os
import ycm_core


def FlagsForFile( filename ):
  dirname = os.path.dirname(filename)
  flags = {'flags': ['--std=c++11', '-MMD', '-Wall', '-I', '.','-I','..','-I'+dirname,'/fwdpp'], 
      'do_cache': True}
  return flags
