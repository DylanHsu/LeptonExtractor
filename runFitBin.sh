outputDir=$1
histName=$2
dataFileName=$3
templateFileName=$4
signalModel=$5
bkgModel=$6
plotTitle=$7
origDir=`pwd`

#run
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home/dhsu/CMSSW_8_0_26_patch1/src
eval `scramv1 runtime -sh`
cd LeptonExtractor
mkdir -p $outputDir
root -b -l <<EOF
assert(0==gSystem->Load("fitBin_C.so"))
printf("fitBin(\"$histName\",\"$dataFileName\",\"$templateFileName\",$signalModel,$bkgModel,\"$plotTitle\",\"$outputDir\")\n");
fitBin("$histName","$dataFileName","$templateFileName",$signalModel,$bkgModel,"$plotTitle","$outputDir")
.q
EOF
cd $origDir
