import os, numpy as np
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-M", "--mtot-inj", dest="M_inj", help='total injected mass')
(options, args) = parser.parse_args()

M_inj = options.M_inj

if M_inj == None:
	print 'cool'


