# PFCalEE

Geant4 simulation of a Si-base sampling calorimeter

Check https://twiki.cern.ch/twiki/bin/view/CMS/HGCalPerformanceStudiesWithG4

Geometry implementation is instantiated by detector versions in an enum - cf. src/DetectorConstruction.cc and src/SamplingSection.cc

A small ntuple is stored with the energy deposits - cf. src/EventAction.cc 

When changing the ttree content, adding homemade classes rather than
simple objects, the classes should be added to userlib. The dictionary
for root to understand the classes also need to be remade. Use "make
dictionary" before make, inside of userlib/. Follow instructions in
https://twiki.cern.ch/twiki/bin/view/Sandbox/AnnemarieMagnanSandbox
for what needs to be put in the homemade classes and makefile for root to
understand them (see also example in class userlib/include/HGCSSSimHit.hh).

## Setup the environment (SLC6)
tcsh
source g4env4lpc.csh
./compile.sh


## For loop to run several jobs in parallel.  Set p = 0, f = 1 for first pass fast events (only events w/ hard nuclear interactions are generated).  Set p = 1 second iteration to select out the seeded events:
# Setting m = 2, v = 5 should give tracker+ecal+hcal -- needs to be verified
for i in `seq 0 5`; do python submitProdLPC_v2.py  -r ${i} -o /afs/cern.ch/work/o/ocolegro/test/ -n 1000 -p 0 -f 1 -m 2 - v 5; done

##Set p = 1 second iteration to select out the seeded events:

for i in `seq 0 5`; do python submitProdLPC_v2.py  -r ${i} -o /afs/cern.ch/work/o/ocolegro/test/ -n 1000 -p 1 -f 1 -m 2 - v 5; done

##To generate full events in one pass set p = 0 and f = 0

for i in `seq 0 5`; do python submitProdLPC_v2.py  -r ${i} -o /afs/cern.ch/work/o/ocolegro/test/ -n 1000 -p 0 -f 0 -m 2 - v 5; done


## For loop is to generate several samples with same stat in parallel.  Set p = 0, f = 1 for first pass fast events (only events w/ hard nuclear interactions are generated).  Set p = 1 second iteration to select out the seeded events:



## use option -S to test the script locally
## to write to eos, one might replace -o /afs/cern.ch/work/o/ocolegro/test/ with -e /store/user/ocolegro/test.
