#!/usr/bin/env python

import os
import sys
import shutil
from argparse import ArgumentParser
import ROOT
import subprocess
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libPandaTreeObjects.so')
cmsswBase = os.environ['CMSSW_BASE'];
os.chdir(cmsswBase+'/src/LeptonExtractor')
ROOT.gROOT.LoadMacro(cmsswBase+'/src/LeptonExtractor/pandaSelection.C+')

CATALOGDIR = '/home/cmsprod/catalog/t2mit/pandaf/005'
TASKNAME = 'tpskim'

argParser = ArgumentParser(description = "Parameters for Anaysis")
argParser.add_argument('task', metavar = 'task',  help = 'Task name.')
argParser.add_argument('dataset', metavar = 'dataset',  help = 'Dataset Name')
argParser.add_argument('fileset', metavar = 'fileset', nargs = '?', help = 'Fileset Name')
argParser.add_argument('--do_electrons', '-E', action = 'store_true', dest = 'do_electrons', help = 'Call for processing of electrons')
argParser.add_argument('--do_muons', '-M', action = 'store_true', dest = 'do_muons', help = 'Call for processing of muons')
argParser.add_argument('--real_data', '-RD', action = 'store_true', dest = 'real_data', help = 'DATA-type files to be processed')
argParser.add_argument('--verbose', '-V', action = 'store_true', dest = 'verbose', help = 'Call for DEBUG level printouts')
argParser.add_argument('--truth_matching', '-TM', action = 'store_true', dest = 'truth_matching', help = 'Call for Monte Carlo Truth Matching')
argParser.add_argument('--max_entries', '-ME', metavar = 'ME', dest = 'max_entries', type = int, default = -1, help = 'Restrict the number of processed events')
argParser.add_argument('--electron_trigger', '-ET', metavar = 'ET', dest = 'electron_trigger', type = int, default = 3, help = 'Electron Trigger value')
argParser.add_argument('--muon_trigger', '-MT', metavar = 'MT', dest = 'muon_trigger', type = int, default = 6, help = 'Muon Trigger value')
argParser.add_argument('--truth_matching_dR', '-dR', metavar = 'dR', dest = 'truth_matching_dR', type = float, default = 0.3, help = 'Truth Matching Delta-R value')
argParser.add_argument('--tag_id', '-TID', metavar = 'TID', dest = 'tag_id', type = int, default = 6, help = 'Tag_Id value')
argParser.add_argument('--tag_iso', '-TIS', metavar = 'TIS', dest = 'tag_iso', type = int, default = 6, help = 'Tag_Iso value')
argParser.add_argument('--probe_id', '-PID', metavar = 'PID', dest = 'probe_id', type = int, default = 0, help = 'Probe_Id value')
argParser.add_argument('--probe_iso', '-PIS', metavar = 'PIS', dest = 'probe_iso', type = int, default = 0, help = 'Probe_Iso value')
argParser.add_argument('--passing_probe_id', '-PPID', metavar = 'PPID', dest = 'passing_probe_id', type = int, default = 6, help = 'Passing_Probe_Id value')
argParser.add_argument('--passing_probe_iso', '-PPIS', metavar = 'PPIS', dest = 'passing_probe_iso', type = int, default = 6, help = 'Passing_Probe_Iso value')
argParser.add_argument('--tag_eta_max', '-TEM', metavar = 'TEM', dest = 'tag_eta_min', type = float, default = 2.1, help = 'Maximum Eta value for tag')
argParser.add_argument('--tag_pt_min', '-TPM', metavar = 'TPM', dest = 'tag_pt_min', type = float, default = 30, help = 'Minimum Pt value for tag')
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

    #outdir = '/mnt/hadoop/scratch/' + os.environ['USER'] + '/' + TASKNAME + '/' + args.dataset + str(args.tag_id)
    outdir = '/mnt/hadoop/scratch/' + os.environ['USER'] + '/' + TASKNAME + '/' + args.dataset
    try:
        os.makedirs(outdir)
    except OSError:
        pass

    submitter.pre_args = 'skim ' + args.dataset
    submitter.post_args = ' '.join(sys.argv[3:])

    filesets = []
    with open(CATALOGDIR + '/' + args.dataset + '/Filesets') as catalog:
    #with open(CATALOGDIR + '/' + args.dataset + '/RawFiles.00') as catalog:
        for line in catalog:
            filesets.append(line.split()[0])
            # break

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
    
    tmpdir = '/tmp/' + os.environ['USER'] + '/' + TASKNAME + '/' + args.dataset
    try:
        os.makedirs(tmpdir)
    except OSError:
        pass

    #outdir = '/mnt/hadoop/scratch/' + os.environ['USER'] + '/' + TASKNAME + '/' + args.dataset + str(args.tag_id)
    outdir = '/mnt/hadoop/scratch/' + os.environ['USER'] + '/' + TASKNAME + '/' + args.dataset

    tmpnameE = tmpdir + '/'  + args.fileset + 'E.root'
    tmpnameM = tmpdir + '/' + args.fileset + 'M.root'
    finalnameE = outdir + '/' + args.fileset + 'E'+'_tag'+str(args.tag_id)+'-'+str(args.tag_iso)+'_probe'+str(args.passing_probe_id)+'-'+str(args.passing_probe_iso)+'.root'
    finalnameM = outdir + '/' + args.fileset + 'M'+'_tag'+str(args.tag_id)+'-'+str(args.tag_iso)+'_probe'+str(args.passing_probe_id)+'-'+str(args.passing_probe_iso)+'.root'

    ROOT.gettagprobe.tag_id = args.tag_id
    ROOT.gettagprobe.tag_iso = args.tag_iso
    ROOT.gettagprobe.probe_id = args.probe_id
    ROOT.gettagprobe.probe_iso = args.probe_iso
    ROOT.gettagprobe.passing_probe_id = args.passing_probe_id
    ROOT.gettagprobe.passing_probe_iso = args.passing_probe_iso
    ROOT.gettagprobe.tag_eta_min = args.tag_eta_min
    ROOT.gettagprobe.tag_pt_min = args.tag_pt_min
  
    

    file_listE = ""
    file_listM = ""
    files = []
    for i in range(len(paths)):
        
        filename = "tnp_" + str(paths[i].rsplit("/",1)[1].replace(".root",""))
        ROOT.make_tnp_skim(paths[i], filename, args.do_electrons, args.do_muons, args.real_data, args.verbose, args.truth_matching, args.max_entries, args.electron_trigger, args.muon_trigger, args.truth_matching_dR)
        files.append(filename)
    
    for name in files:
        
        file_listE = file_listE + os.getcwd() + "/%s " % (name + "_electronTnP.root")
        file_listM = file_listM + os.getcwd() + "/%s " % (name + "_muonTnP.root")

    if args.do_electrons:
        subprocess.call("cd %s/src; cmsenv; hadd -f %s %s"  % (os.environ['CMSSW_BASE'],tmpnameE, file_listE), shell = True)
        shutil.copy(tmpnameE, finalnameE)
        os.remove(tmpnameE)

    if args.do_muons:    
        subprocess.call("cd %s/src; cmsenv; hadd -f %s %s"  % (os.environ['CMSSW_BASE'],tmpnameM, file_listM), shell = True)
        shutil.copy(tmpnameM, finalnameM)
        os.remove(tmpnameM)
