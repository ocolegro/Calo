#!/usr/bin/env python

import os
absPath = '/uscms/home/ocolegro/Calo/signal/4GeV_nodecay_modified/'
outDir  = '/store/user/ocolegro/4GeV_nodecay_magOn_run2/'
files = os.listdir('signal/4GeV_nodecay_modified')
#python submitProdLPC_signal.py -f /uscms/home/ocolegro/Calo/signal/4GeV/SLAC.4GeV.W.mchi.0.05.map.0.15.alpha0.1.fermionDM_unweighted_events_modified.lhe -s 1 -e /eos/uscms/store/user/ocolegro/signal
for file in files:
    	print 'Calling %s ' % ('python submitProdLPC_signal.py -f %s/%s -s 0 -e %s' % (absPath,file,outDir))
	os.system('python submitProdLPC_signal.py -f %s/%s -s 0 -e %s -n 100000' % (absPath,file,outDir))
