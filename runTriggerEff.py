#!/usr/bin/env python

import os
import sys
import shutil
from argparse import ArgumentParser
import ROOT
import subprocess
ROOT.gROOT.SetBatch(True)
#ROOT.gSystem.Load('libPandaTreeObjects.so')
cmsswBase = os.environ['CMSSW_BASE'];
#os.chdir(cmsswBase+'/src/LeptonExtractor')
ROOT.gROOT.LoadMacro(cmsswBase+'/src/LeptonExtractor/triggerEff.C+')

#CATALOGDIR = '/home/cmsprod/catalog/t2mit/pandaf/009'
CATALOGDIR = '/home/cmsprod/catalog/t2mit/pandaf/010'
TASKNAME = 'triggerEff'

# LeptonExtractor/runTriggerEff.py submit SingleElectron+Run2017B-31Mar2018-v1+MINIAOD SETrigSoup --real_data

argParser = ArgumentParser(description = "Parameters for Anaysis")
argParser.add_argument('task', metavar = 'task',  help = 'Task name.')
argParser.add_argument('dataset', metavar = 'dataset',  help = 'Dataset Name')
argParser.add_argument('trigtype', metavar = 'dataset',  help = 'Trigger Type')
argParser.add_argument('fileset', metavar = 'fileset', nargs = '?', help = 'Fileset Name')
argParser.add_argument('--real_data', '-RD', action = 'store_true', dest = 'real_data', help = 'DATA-type files to be processed')
argParser.add_argument('--debug', '-v', action = 'store_true', dest = 'debug', help = 'Debug mode')
args = argParser.parse_args()

if args.task == 'submit':
    sys.path.append(cmsswBase+'/src/LeptonExtractor/condor_run_lib')
    from condor_run import CondorRun
    
    submitter = CondorRun(os.path.realpath(__file__))
    submitter.logdir = '/local/' + os.environ['USER']
    submitter.hold_on_fail = True
    submitter.min_memory = 1

    if not os.path.exists(submitter.logdir):
        os.makedirs(submitter.logdir)

    outdir = '/mnt/hadoop/scratch/' + os.environ['USER'] + '/' +TASKNAME+'/'+args.trigtype+'/'+args.dataset
    try:
        os.makedirs(outdir)
    except OSError:
        pass

    submitter.pre_args = 'skim ' + args.dataset + ' ' + args.trigtype
    submitter.post_args = ' '.join(sys.argv[4:])

    filesets = []
    with open(CATALOGDIR + '/' + args.dataset + '/Filesets') as catalog:
        for line in catalog:
            filesets.append(line.split()[0])

    submitter.job_args = filesets
        
    submitter.submit(name = TASKNAME)

elif args.task == 'skim':
    with open(CATALOGDIR + '/' + args.dataset + '/Filesets') as catalog:
        for line in catalog:
            if line.startswith(args.fileset):
                #directory = line.split()[1].replace('root://xrootd.cmsaf.mit.edu/', '/mnt/hadoop/cms')
                directory = line.split()[1]
                #break

    paths = []

    with open(CATALOGDIR + '/' + args.dataset + '/Files') as catalog:
        for line in catalog:
            if line.startswith(args.fileset):
                paths.append(directory + '/' + line.split()[1])
    
    tmpdir = '/tmp/' + os.environ['USER']+'/'+TASKNAME+'/'+args.trigtype+'/'+args.dataset
    try:
        os.makedirs(tmpdir)
    except OSError:
        pass

    outdir = '/mnt/hadoop/scratch/'+os.environ['USER']+'/'+TASKNAME+'/'+args.trigtype+'/'+args.dataset

    tmpname = tmpdir + '/'  + args.fileset + '.root'
    finalname = outdir + '/' + args.fileset + '.root'

    file_list = ""
    files = []
    for i in range(len(paths)):
        filename = "triggerEff_" + str(paths[i].rsplit("/",1)[1])
        ROOT.triggerEff(args.trigtype, paths[i], filename, args.real_data, args.debug) 
        files.append(filename)
    
    for name in files:
        file_list = file_list + os.getcwd() + "/%s " % (name)

    subprocess.call("cd %s/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval `scramv1 runtime -sh`; hadd -f %s %s"  % (os.environ['CMSSW_BASE'],tmpname, file_list), shell = True)
    shutil.copy(tmpname, finalname)
    os.remove(tmpname)

