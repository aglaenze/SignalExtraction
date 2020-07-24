#!/bin/bash
clean()
{
rm *.so
rm *.d
rm *.pcm
}

clean

muonfilter=std
#muonfilter=nopXdca
#muonfilter=noLpt

#I have a diffrent folder inside "rootFiles/" for each of the muon filter processed in the analysis
#in alidock
path_to_rootfiles_data=\"/home/alidock/analysis-alice/p-Pb-2016/rootFiles/$muonfilter\"
path_to_rootfiles_MC=\"/home/alidock/analysis-alice/p-Pb-2016/rootFiles\"
#outside alidock
#path_to_rootfiles_data=\"/Users/aglaenzer/alidock/analysis-alice/p-Pb-2016/rootFiles/$muonfilter\"
#path_to_rootfiles_MC=\"/Users/aglaenzer/alidock/analysis-alice/p-Pb-2016/rootFiles\"

# Also, I asume that the format of the root file names is
# - AnalysisResults_LHC16r_MC_kIncohJpsiToMu.root for MC data
# - AnalysisResults_LHC16r.root for real data

periods="{\"LHC16r\", \"LHC16s\"}"
#periods="{\"LHC16r\"}"


mMin=2.5
mMax=3.9
ptMin=0
ptMax=3.5
#useCuts=true
useCuts=true
logScale=true
drawPulls=false		# draws graphs data-fit, in any case the graphs (data-fit)/sigma are plotted


echo $muonfilter
echo $periods
root -l -q "ToyMC.C+($path_to_rootfiles_MC, $drawPulls)"
clean
root -l -q "TailParameters.C+($path_to_rootfiles_MC, $periods, $logScale)"
clean
root -l -q "Splot.C+($path_to_rootfiles_data, $periods, $mMin, $mMax, $ptMin, $ptMax, $useCuts)"
clean
root -l -q "TwoDPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods, $mMin, $mMax, $ptMin, $ptMax, $useCuts, $logScale, $drawPulls)"
clean
