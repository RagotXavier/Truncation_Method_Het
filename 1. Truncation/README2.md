The code solves the steady state in Julia. Then, it uses the DYNARE solver for the dynamics (perturbation around the steady state). DYNARE is also used to double-check the steady state. DYNARE needs MATLAB to run.

The steady state of the model uses the Endogenous Grid Method (EGM), using the code of Alidstair McKay (All errors are ours) 


I - Julia code to solve for the steady state :**

- Main file « Main\_Truncation.jl » It uses the following files:     
  - « Aiyagari\_solve.jl » solve aiyagari model
  - « Parameters.jl » parameters of the Aiyagari model
  - « Projection\_Truncation.jl » truncating the model
  - ` `Save relevant variables in « todynare\_Truncation.mat »

II - DYNARE simulations :** 

- « Dynare\_Truncation.m », write the DYNARE file for full set of instruments and launch it. Save results in
  - « tofigtruncation.mat » impulse response functions from the variables 
  - « tofigtruncationsim.mat » time series simulations for the main variables 
  - « tofigtruncationsimt.mat » time series simulations for the main variables (deviation from the steady state)
  - « tofigtruncation.mat » results of the dynamics


