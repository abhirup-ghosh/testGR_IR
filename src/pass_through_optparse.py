#!/usr/bin/python

# "Pass-through" option parsing -- an OptionParser that ignores
# unknown options and lets them pile up in the leftover argument
# list.  Useful for programs that pass unknown options through
# to a sub-program.

from optparse import OptionParser, BadOptionError

class PassThroughOptionParser(OptionParser):

    def _process_long_opt(self, rargs, values):
        try:
            OptionParser._process_long_opt(self, rargs, values)
        except BadOptionError, err:
            self.largs.append(err.opt_str)

    def _process_short_opts(self, rargs, values):
        try:
            OptionParser._process_short_opts(self, rargs, values)
        except BadOptionError, err:
            self.largs.append(err.opt_str)
