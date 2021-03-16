"""
Functions to convert uvfits to ms with CASA

Author: Jack Line <https://github.com/JLBLine>
"""

from subprocess import call
from sys import argv

import numpy as np


def make_ms(uvfits_file, obs):
    name = uvfits_file[:-7]

    call("rm -r %s.ms" % name, shell=True)
    importuvfits(fitsfile=uvfits_file, vis="%s.ms" % name)

    # call('fixmwams %s.ms %s' %(name,obs),shell=True)

    # print('/usr/local/cotter/build/fixmwams %s.ms metafits/%s_metafits_ppds.fits' %(name,obs))


def add_MWA_name(name):
    msfile = "%s.ms" % name
    tb.open(msfile + "/OBSERVATION", nomodify=False)
    tb.putcol("TELESCOPE_NAME", "MWA")
    tb.flush()


for band in np.arange(1, 25):
    make_ms("%s%02d.uvfits" % (argv[3], band), argv[4])
    add_MWA_name("%s%02d" % (argv[3], band))
