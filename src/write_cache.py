import os, re, commands

def write_cache(frame, cache_folder=None, cachefilename=None, channel=None, datalen=None, formatter=':'):
  """
  Writes cache for frame with a single channel.
  
  Parameters
  ----------
  frame : frame file
  cache_folder : output folder to write file to
  
  Returns
  -------
  A blank-separated string (needs to be split):
  'cachefilename ch channel start_time datalen realpath'
  """
  frame_dump = commands.getoutput('FrDump -i %s'%(frame))
  start_time = int(re.findall('start at:[0-9]* ', frame_dump)[0][9:-1])
  if datalen is None:
    datalen = int(re.findall('length=[0-9]*\.', frame_dump)[0][7:-1])
  if channel is None:
    channel = re.findall('ProcData:.*', frame_dump)[0][10:]
  realpath = os.path.realpath(frame)
  if cachefilename is None:
    cachefilename = '%s-%s_CACHE-%d-%d.lcf'%(channel[0], channel[3:], start_time, datalen+16)
  if cache_folder is not None:
    os.system('mkdir -p %s'%(cache_folder))
    ofile = open(os.path.join(cache_folder, cachefilename), 'w')
    ofile.write('%s %s%s%s %d %d file://localhost%s\n'%(channel[0], channel[:2], formatter, channel[3:], start_time, datalen, realpath))
    ofile.close()
  return '%s %s %s %d %d file://localhost%s'%(cachefilename, channel[0], channel, start_time, datalen, realpath)

