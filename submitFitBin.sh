jobDir=$CMSSW_BASE/src/LeptonExtractor/$1
cacheDir=$2
cd $CMSSW_BASE/src/LeptonExtractor
mkdir -p $jobDir/plots 2>/dev/null
#rm fitBin_C.d  fitBin_C.so  fitBin_C_ACLiC_dict_rdict.pcm 2>/dev/null
#root -b -l <<EOF
#.L fitBin.C+g
#EOF
cd ..
if [ "$cacheDir" == "" ]
then
  #cacheDir="$jobDir/submit"
  #mkdir -p $cacheDir
  PandaCore/bin/submit --exec LeptonExtractor/runFitBin.sh --arglist $jobDir/jobArgs.txt
else
  PandaCore/bin/submit --exec LeptonExtractor/runFitBin.sh --arglist $jobDir/jobArgs.txt --cache $cacheDir
fi

htmlOut=$jobDir/plots.html
rm $jobDir/*.html 2>/dev/null
echo "<html> " >> $htmlOut
echo "<body bgcolor=\"EEEEEE\"><h1>Efficiency extraction fits</h1>" >> $htmlOut
echo "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" >> $htmlOut
i=0
while read -r line
do
  if [ "$((i % 4))" == "0" ]; then echo "<tr>" >> $htmlOut; fi
  words=( $line )
  echo "<td width=\"25%\"><a target=\"_blank\" href=\"plots/fitres_${words[1]}.txt\"><img src=\"plots/${words[1]}.png\" alt=\"plots/${words[1]}.png\" width=\"100%\"></a></td>" >> $htmlOut
  if [ "$((i % 4))" == "3" ]; then echo "</tr>" >> $htmlOut; fi
  i=$((i+1))
done < "$jobDir/jobArgs.txt"
if [ "$((i % 4))" != "0" ]; then echo "<td width=\"25%\"></td><td width=\"25%\"></td></tr>" >> $htmlOut; fi
echo "</table>" >> $htmlOut
echo "<hr />" >> $htmlOut
echo "</body>" >> $htmlOut
echo "</html>" >> $htmlOut
#  unsigned int i;
#  for(i=0; i<nbins; i++) {
#    if(i%2==0) echo "<tr>"     >> $htmlOut
#    echo "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pass" << name << "_" << i << ".png\"><img src=\"plots/pass" << name << "_" << i << ".png\"alt=\"plots/pass" << name << "_" << i << ".png\" width=\"100%\"></a></td>" >> $htmlOut
#    echo "<td width=\"25%\"><a target=\"_blank\" href=\"plots/fail" << name << "_" << i << ".png\"><img src=\"plots/fail" << name << "_" << i << ".png\"alt=\"plots/fail" << name << "_" << i << ".png\" width=\"100%\"></a></td>" >> $htmlOut
#    if(i%2) echo "</tr>" >> $htmlOut
#  }
#  if(i%2) {
#    echo "<td width=\"25%\"></td>" >> $htmlOut
#    echo "<td width=\"25%\"></td>" >> $htmlOut
#    echo "</tr>" >> $htmlOut
#  }
#  echo "</table>" >> $htmlOut
#  echo "<hr />" >> $htmlOut
#    
#  echo "</body>" >> $htmlOut
#  echo "</html>" >> $htmlOut

