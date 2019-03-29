"""
Part of header common to all projects.

(C) Archisman Ghosh, 2014-Oct-29
"""
import commands

machine = commands.getoutput('hostname')
user = commands.getoutput('whoami')
home = commands.getoutput('echo ${HOME}')
if user=='ajith':
  work_dir = '%s/working'%(home)
elif user=='archis' or user=='carchisman':
  work_dir = '%s/Work'%(home)
elif user=='wdp':
  work_dir = '%s/Dropbox'%(home)
elif user=='wdpozzo':
  work_dir = '%s'%(home)
elif user=='nishad' or user=='nishadb':
  work_dir = '%s/GW_Projects'%(home)
elif user=='sudarshan':
  work_dir = '%s/Documents/ICTS/Work/pyPE'%(home)
