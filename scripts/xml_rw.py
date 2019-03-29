#!/usr/bin/env python

"""
Script to modify specific columns of sim_inspiral:table in xml file created by lalapps_inspinj. Implemented here only for z component of spin.

(C) Archisman Ghosh, 2015-10-23
"""

import os, re, optparse

parser = optparse.OptionParser()
parser.add_option('-i', '--infile', type='string', help='input file name')
parser.add_option('-o', '--outfile', type='string', default=None, help='output file name')
parser.add_option('--spin1x', type='float', help='spin1x to write')
parser.add_option('--spin1y', type='float', help='spin1y to write')
parser.add_option('--spin1z', type='float', help='spin1z to write')
parser.add_option('--spin2x', type='float', help='spin2x to write')
parser.add_option('--spin2y', type='float', help='spin2y to write')
parser.add_option('--spin2z', type='float', help='spin2z to write')
(options, args) = parser.parse_args()

spin1x = options.spin1x
spin1y = options.spin1y
spin1z = options.spin1z
spin2x = options.spin2x
spin2y = options.spin2y
spin2z = options.spin2z
infilename = options.infile
outfilename = options.outfile
overwrite = False

# Set a temporary output filename if one is not provided
if outfilename is None or outfilename==infilename:
  outfilename = infilename[:-4]+'_mod.xml'
  overwrite = True

# Set some flags that tell you which part of the xml file are you in
INSIDE_SIM_INSPIRAL_HEADER = False
INSIDE_SIM_INSPIRAL_STREAM = False
colnum = 0

# Open the stream to write in
ofile = open(outfilename, 'w')

# Scan though the input file line by line
for (line_num, line_text) in enumerate(open(infilename, 'r')):
  
  # Enter sim_inspiral:table header
  if re.findall('<Table Name="sim_inspiralgroup:sim_inspiral:table"', line_text):
    INSIDE_SIM_INSPIRAL_HEADER = True
  
  # Exit sim_inspiral:table header and enter sim_inspiral:table stream
  elif re.findall('<Stream Name="sim_inspiralgroup:sim_inspiral:table"', line_text):
    INSIDE_SIM_INSPIRAL_HEADER = False
    INSIDE_SIM_INSPIRAL_STREAM = True
  
  # Scan through the column names in the header to find out column numbers to be changed
  # NOTE: Add more columns here if necessary
  elif INSIDE_SIM_INSPIRAL_HEADER:
    colname = re.findall('Name="sim_inspiralgroup:sim_inspiral:.*?"', line_text)[0][37:-1]
    #print 'column %d:\t%s'%(colnum, colname)
    if colname == 'spin1z':
      spin1z_colnum = colnum
    elif colname == 'spin2z':
      spin2z_colnum = colnum
    elif colname == 'spin1x':
      spin1x_colnum = colnum
    elif colname == 'spin2x':
      spin2x_colnum = colnum
    elif colname == 'spin1y':
      spin1y_colnum = colnum
    elif colname == 'spin2y':
      spin2y_colnum = colnum
    colnum += 1 
  
  # Change the appropriate columns in the stream
  elif INSIDE_SIM_INSPIRAL_STREAM:
    
    # Exit sim_inspiral:table stream if end is hit
    if re.findall('</Stream>', line_text):
      INSIDE_SIM_INSPIRAL_STREAM = False
    
    # Make the actual changes: split, change, join back
    # NOTE: Change more columns here if necessary
    else:
      line_text_split = line_text.split(',')
      line_text_split[spin1x_colnum] = '%e'%(spin1x)
      line_text_split[spin2x_colnum] = '%e'%(spin2x)
      line_text_split[spin1y_colnum] = '%e'%(spin1y)
      line_text_split[spin2y_colnum] = '%e'%(spin2y)
      line_text_split[spin1z_colnum] = '%e'%(spin1z)
      line_text_split[spin2x_colnum] = '%e'%(spin2x)
      line_text_split[spin2y_colnum] = '%e'%(spin2y)
      line_text_split[spin2z_colnum] = '%e'%(spin2z)
      line_text = ','.join(line_text_split)
  
  # Write the line back to the output stream
  ofile.write(line_text)

# Close the stream
ofile.close()

# Overwrite the new file to the old file
if overwrite:
  os.system('mv %s %s'%(outfilename, infilename))
