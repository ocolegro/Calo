#!/usr/bin/env python

import os
absPath = '/uscms/home/ocolegro/Calo/signal/4GeV/'
outDir  = '/eos/uscms/store/user/ocolegro/signalInc_test'
files = os.listdir('signal/4GeV')
#python submitProdLPC_signal.py -f /uscms/home/ocolegro/Calo/signal/4GeV/SLAC.4GeV.W.mchi.0.05.map.0.15.alpha0.1.fermionDM_unweighted_events_modified.lhe -s 1 -e /eos/uscms/store/user/ocolegro/signal
for file in files[0:2]:
    os.system('python submitProdLPC_signal.py -f %s/%s -s 0 -e %s' % (absPath,file,outDir))
