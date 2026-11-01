# Plots for ekman-astz, kept out of the compositing scripts so their worker processes
# never load PythonPlot/CondaPkg. Runs in the plot/ environment, from plot/:
#   julia --project=. plot_bsiso_phase.jl
using DelimitedFiles
using PythonPlot

# BSISO1 and BSISO2 phase time series, from BSISO_phase.txt written by bsiso_phase.jl
# YEAR DOY MONTH DATE BSISO1_phase BSISO2_phase
P = readdlm("../BSISO_phase.txt", Int, skipstart=1)
plot( P[:,5] ) # increases
plot( P[:,6] )
