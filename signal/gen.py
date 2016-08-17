#!/usr/bin/env python

import os
files = os.listdir('4GeVInclusive')
for file in files:
    os.system('python LHEParser.py -f 4GeVInclusive/%s -x 0 -y 0 -z -25.75' % (file))
