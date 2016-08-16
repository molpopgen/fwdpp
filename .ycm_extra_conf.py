import os
import ycm_core


def FlagsForFile( filename ):
  dirname = os.path.dirname(filename)
  flags = {'flags': ['--std=c++11', '-MMD', '-Wall', '-I', '.', '-I' '../include', '-I', dirname, '-I', dirname +
      '/../include',
      '-I',dirname+'/../../include'], 'do_cache': True}
  return flags
