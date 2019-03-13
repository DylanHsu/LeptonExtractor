outputDir=$1
histName=$2
dataFileName=$3
sigTemplateFileName=$4
bkgTemplateFileName=$5
bkgTemplateFileName2=$6
signalModel=$7
bkgModel=$8
plotTitle=$9
signalLabel=${10}
bkgLabel=${11}
origDir=`pwd`

#run
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home/dhsu/CMSSW_9_4_6/src
eval `scramv1 runtime -sh`
cd LeptonExtractor
mkdir -p $outputDir
#assert(0==gSystem->Load("fitBin_C.so"))
root -b -l <<EOF
RooFitResult *fr=0;
gSystem->Load("fitBin_C.so")
printf("fitBin(\"$histName\",\"$dataFileName\",\"$sigTemplateFileName\",\"$bkgTemplateFileName\",\"$bkgTemplateFileName2\",$signalModel,$bkgModel,\"$plotTitle\",\"$outputDir\",\"$signalLabel\",\"$bkgLabel\")\n");
fr=fitBin("$histName","$dataFileName","$sigTemplateFileName","$bkgTemplateFileName","$bkgTemplateFileName2",$signalModel,$bkgModel,"$plotTitle","$outputDir","$signalLabel","$bkgLabel");
assert(fr);
.q
EOF
cd $origDir
