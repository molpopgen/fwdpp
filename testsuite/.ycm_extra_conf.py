import os
import ycm_core


def FlagsForFile( filename ):
  dirname = os.path.dirname(filename)
  flags = {'flags': ['--std=c++11', '-Wall', '-I', '.','-I','..','-I','../..'], 
      'do_cache': True}
  return flags
