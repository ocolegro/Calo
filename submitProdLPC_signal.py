#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
#Long Que is 2nw
parser.add_option('-q', '--queue' ,    dest='queue'             , help='batch queue'            , default='1nd')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='hexaV02-01-01')
parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,      type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=2,      type=int)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=1,      type=int)
parser.add_option('-g', '--generator'      ,    dest='generator'             , help='generator flag'                  , default=1,      type=int)
parser.add_option('-s', '--speed'      ,    dest='speed'             , help='speed flag'                  , default=0,      type=int)
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to generate' , default=10000,    type=int)
parser.add_option('-f', '--file'       ,    dest='file'              , help='file name' , default='')
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

label=''


nevents=opt.nevts
myqueue=opt.queue

print 'creating the job'
#eosDir='%s/%s'%('','mchi_' +fname.split('.')[4] + "." +  fname.split('.')[5] + '_GEV_alpha_' + fname.split('.')[7] + "." + fname.split('.')[8]   + '_GEV')

#eosDir='%s/%s'%(opt.eos,'mchi_' +opt.file.split('.')[4] + "." +  opt.file.split('.')[5] + '_GEV_alpha_' + opt.file.split('.')[7] + "." + opt.file.split('.')[8]   + '_GEV')
#massJob = (opt.file.split('//')[-1].split('_')[0].split('.map')[-1]).split('.')
massJob = (opt.file.split('//')[-1].split('_')[0].split('.map')[-1])

eosDir='%s/%s'%(opt.eos,'M_alpha_%sGEV_%s' % ( massJob[0:len(massJob)-2], massJob[-1]))

print 'xrdfs root://cmseos.fnal.gov mkdir %s'%eosDir
#os.system('xrdfs root://cmseos.fnal.gov rm %s'%eosDir)
os.system('xrdfs root://cmseos.fnal.gov mkdir %s'%eosDir)
os.system('xrdfs root://cmseos.fnal.gov rm  /%s/PFCalEE' % eosDir)
os.system('xrdfs root://cmseos.fnal.gov rm  /%s/g4env4lpc.sh' % eosDir)
os.system('xrdfs root://cmseos.fnal.gov rm  /%s/libPFCalEE.so' % eosDir)
os.system('xrdfs root://cmseos.fnal.gov rm  /%s/libPFCalEEuserlib.so' % eosDir)
os.system('xrdfs root://cmseos.fnal.gov rm  /%s/runJob.sh' % eosDir)
os.system('xrdfs root://cmseos.fnal.gov rm  /%s/g4steer.mac' % eosDir)
os.system('xrdfs root://cmseos.fnal.gov rm  /%s/submit.jdl' % eosDir)
os.system('xrdfs root://cmseos.fnal.gov rm  /%s/b18d36.dat' % eosDir)
os.system('xrdfs root://cmseos.fnal.gov rm  %s/%s' % (eosDir,opt.file))


#print 'The outdir is %s' % eosDir
os.system('xrdcp $HOME/geant4_workdir/bin/Linux-g++/PFCalEE root://cmseos.fnal.gov/%s/' % eosDir)
#os.system('xrdcp g4env4lpc.sh root://cmseos.fnal.gov/%s/' % eosDir)
os.system('xrdcp b18d36.dat root://cmseos.fnal.gov/%s/' % eosDir)
os.system('xrdcp $HOME/geant4_workdir/tmp/Linux-g++/PFCalEE/libPFCalEE.so root://cmseos.fnal.gov/%s/' % eosDir)
os.system('xrdcp userlib/lib/libPFCalEEuserlib.so root://cmseos.fnal.gov/%s/' % eosDir)
os.system('xrdcp %s root://cmseos.fnal.gov/%s/' % (opt.file,eosDir))


#wrapper
scriptFile = open('runJob.sh', 'w')
scriptFile.write('#!/bin/bash\n')
scriptFile.write('source g4env4lpc.sh\n')#%(os.getcwd()))

scriptFile.write('./PFCalEE g4steer.mac %d %d %d %d | tee g4.log\n'%(opt.version,opt.model,opt.generator,opt.speed))

scriptFile.write('xrdcp -f PFcal.root root://cmseos.fnal.gov/%s/%s.root\n'%(eosDir,opt.file.split('/')[-1].replace('.lhe','')))
scriptFile.write('xrdfs root://cmseos.fnal.gov rm  /%s/PFcal.root\n' % eosDir)

scriptFile.write('localdir=`pwd`\n')
scriptFile.write('echo "--Local directory is " $localdir >> g4.log\n')
scriptFile.write('ls * >> g4.log\n')



scriptFile.write('echo "--deleting core files: too heavy!!"\n')
scriptFile.write('rm core.*\n')
scriptFile.write('echo "All done"\n')
scriptFile.close()
os.system('xrdcp %s root://cmseos.fnal.gov/%s/' % ('runJob.sh',eosDir))



print 'submitting to the cluster'
#write geant 4 macro
g4Macro = open('g4steer.mac', 'w')
g4Macro.write('/control/verbose 0\n')
g4Macro.write('/control/saveHistory\n')
g4Macro.write('/run/verbose 0\n')
g4Macro.write('/event/verbose 0\n')
g4Macro.write('/tracking/verbose 0\n')
g4Macro.write('/filemode/inputFilename %s\n' % opt.file.split('/')[-1])
g4Macro.write('/random/setSeeds %d %d\n'%( random.uniform(0,100000), random.uniform(0,100000) ) )
g4Macro.write('/run/initialize\n')
g4Macro.write('/run/beamOn %d\n'%(nevents))
g4Macro.close()
os.system('xrdcp %s root://cmseos.fnal.gov/%s/' % ('g4steer.mac',eosDir))

#submit
#os.system('echo %s ' %('chmod 777 %s/runJob.sh'%eosDir))

os.system('chmod 777 //eos/uscms%s/PFCalEE' % (eosDir))
os.system('chmod 777 //eos/uscms%s/runJob.sh' % (eosDir))

if opt.nosubmit : os.system('LSB_JOB_REPORT_MAIL=N echo bsub -q %s -N %s/runJob.sh'%(myqueue,eosDir))
else:
    #os.system("LSB_JOB_REPORT_MAIL=N bsub -q %s -N \'%s/runJob.sh\'"%(myqueue,eosDir))
    name = "submitRun%s" % (opt.run)
    f2n = "submit.jdl" ;
    outtag = "out_%s_$(Cluster)" % (name)
    f2=open(f2n, 'w')
    f2.write("universe = vanilla \n");
    f2.write("Executable = %s \n" % ('//eos/uscms%s/runJob.sh'%(eosDir)) );
    f2.write('Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )\n');
    f2.write("request_disk = 100000\n");
    f2.write("request_memory = 1000\n");
    f2.write("Should_Transfer_Files = YES \n");
    f2.write("Transfer_Input_Files = b18d36.dat,g4env4lpc.sh,libPFCalEE.so,libPFCalEEuserlib.so,PFCalEE,g4steer.mac,%s \n" %(opt.file.split('/')[-1]) );
    f2.write("WhenToTransferOutput  = ON_EXIT_OR_EVICT \n");
    f2.write("Output = "+outtag+".stdout \n");
    f2.write("Error = "+outtag+".stderr \n");
    f2.write("Log = "+outtag+".log \n");
    f2.write("Notification    = Error \n");
    f2.write("x509userproxy = $ENV(X509_USER_PROXY) \n")
    f2.write("Queue 1 \n");
    f2.close();
    os.system('xrdcp %s root://cmseos.fnal.gov/%s/' % ('submit.jdl',eosDir))
    print 'Changing dir to %s' % (eosDir)
    os.chdir("//eos/uscms%s" % (eosDir));
    os.system("condor_submit submit.jdl");# % (submit.jdl));