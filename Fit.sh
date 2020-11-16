#!/bin/bash

## variables
periods="{\"LHC16r\", \"LHC16s\"}"
#periods="{\"LHC16s\"}"

mMin=2.0
mMax=2.8

#mMin=2.7
#mMax=3.4

ptMin=0
ptMax=5
useCuts=1
logScale=0
drawPulls=0			# draws graphs data-fit, in any case the graphs (data-fit)/sigma are plotted
exp=1				# if true, use the exponential for the dissociative contribution. if false, use the power law
exclusiveOnly=0

## end of variables

initiateParameters() {
if test -f $inputfile
then
rm $inputfile
fi

touch $inputfile
echo "
bExc = $bExc
gammaPbYield = $gammaPbYield

mMin = $mMin
mMax = $mMax

ptMin = $ptMin
ptMax = $ptMax

useCuts = $useCuts
logScale = $logScale
drawPulls = $drawPulls		# draws graphs data-fit, in any case the graphs (data-fit)/sigma are plotted
exp = $exp	# if true, use the exponential for the dissociative contribution. if false, use the power law
exclusiveOnly = $exclusiveOnly
" >> $inputfile

}


clean() {
filesToDelete="*.so *.d *.pcm *ACLiC* *.in *.out Include/*.so Include/*.d Include/*.pcm"
for entry in $filesToDelete
do
if test -f $entry;
then
#echo Deleting $entry
rm $entry
fi
done
}


muonfilter=std
#muonfilter=nopXdca
#muonfilter=noLpt

#I have a diffrent folder inside "rootFiles/" for each of the muon filter processed in the analysis
#in alidock
#path="/Users/aglaenzer/alidock/analysis-alice/p-Pb-2016/rootFiles"
#outside alidock
path="/Users/aglaenzer/alidock/analysis-alice/p-Pb-2016/rootFiles"
path_to_rootfiles_data=\"${path}/$muonfilter\"
path_to_rootfiles_MC=\"${path}/MC-$muonfilter\"

# It is assumed that the format of the root file names is
# - AnalysisResults_LHC16r_MC_kIncohJpsiToMu.root for MC data
# - AnalysisResults_LHC16r.root for real data



echo $muonfilter
echo $periods

run() {
clean
root -l -q "Splot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
exit
root -l -q "TwoDPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
#exit
}

run


exp=false

for ((k=0;k<100;k++)); do

# first initiate LHC16r
inputfile="input-LHC16r.txt"
bExc=$(echo "scale=2; 2.8+${k}*0.02" | bc)
bExc=3.5
gammaPbYield=57
initiateParameters

# then LHC16s
inputfile="input-LHC16s.txt"
bExc=$(echo "scale=2; 4.5+${k}*0.02" | bc)
bExc=5.3
gammaPbYield=86
initiateParameters
root -l -q "TwoDPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
done

exit

exclusiveOnly=false
exp=true
run
exit
exp=false
run
exit
exclusiveOnly=false
exp=true
run
exp=false
run
exit



clean
mMin=1.3
mMax=2.3
root -l -q "Background.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods, $mMin, $mMax, $ptMin, $ptMax)"
clean
exit



clean
#root -l -q "GetTemplates.C+($path_to_rootfiles_data, $periods, $mMin, $mMax, $ptMin, $ptMax)"
clean
#exit


root -l -q "NakedPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods, $mMin, $mMax, $ptMin, $ptMax, $useCuts)"
clean
