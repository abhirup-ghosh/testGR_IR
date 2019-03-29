import commands
import numpy as np
#import my_plotter as mp
import matplotlib.pyplot as plt

def load_line(file_name, n=1, **kwargs):
  """
  Reads the first n lines from a file. Default n=1. If n<0, reads the last n lines.
  """
  num_lines = int(commands.getoutput('cat %s | wc -l'%file_name))+int(commands.getoutput('tail -c 1 %s'%file_name)!='')
  if n<0:
    n = num_lines+1+n
  return np.genfromtxt(file_name, skip_header=n-1, skip_footer=num_lines-n, **kwargs)

