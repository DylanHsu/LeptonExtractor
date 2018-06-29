outputDir="/mnt/hadoop/scratch/dhsu/triggerEff/"
trigType=$1
inputFileName=$2
basename=$(basename $2)
outputFileName="${outputDir}/${trigType}/${basename}"
mkdir -p ${outputDir}/${trigType}

origDir=`pwd`
cmsswDir=$CMSSW_BASE/src
tempOutput=$origDir/tempOutput.root
#tempInput=$origDir/tempInput.root
#cp -v $inputFileName $tempInput

#run
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $cmsswDir
mkdir -p $(dirname "${2}")
eval `scramv1 runtime -sh`
root -b -l <<EOF
bool loadedMacro=(0==gSystem->Load("LeptonExtractor/triggerEff_C.so"));
if(!loadedMacro) throw std::runtime_error("Could not load macro shared object, go compile it in ACLiC");
printf("triggerEff(\"$trigType\",\"$inputFileName\",\"$tempOutput\")\n");
triggerEff("$trigType","$inputFileName","$tempOutput");
.q
EOF
cd $origDir
cp -v $tempOutput $outputFileName
rm $tempOutput
