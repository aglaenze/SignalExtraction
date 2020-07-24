# Project for signal extraction using a combined fit in mass and transverse momentum

The simplest is to use ./Fit2d.sh and change inside the folder containing the input root files (for data and Monte Carlo).

- TailParameters.C: used to fit MC mass distribution for J/Psi and Psi(2s), and extract the tail parameters. These values will then be fixed to fit the data.

- SPlot.C: used to get weights on J/Psi, Psi(2s) and background contributions, based on the mass fit. Then it plots pt distribution with weights for J/Psi and Background. The resulting pt distribution for the background and for J/Psi are exported to a rootFile, located at the same place as the input root files (for data).

- TwoDPlot.C: makes the 2D fit of signal within the given mass range and pt range

- ToyMC.C: mixes MC data of J/Psi and $$ \gamma \gamma \rightarrow \mu \mu $$, and then tries sPlot to reconstruct signal. Goal: assess quality of sPlot, given the fact that it should work only for uncorrelated variables. In the case of $$ \gamma \gamma \rightarrow \mu \mu $$, pT and m are correlated, so this aims to test if it fails to reconstruct the 2 contributions.
