This repo contains the main files to implement the truncation method in heterogeneous-agent models. The associated article is “The truncation method to solve heterogeneous-agent models”, by François Le Grand and Xavier Ragot.

1.	First, to understand the algorithm, the best is to start from the file “0. Simple Case/Truncation\_simple\_case.m” This files is self-contained and commented to solve the model with aggregate shocks. It requires Matlab/Octave and Dynare.
 2.	The Folder “1. Truncation” provides additional codes and is more efficient in solving the model. It requires the Julia language 
   a.	In this folder there is a README file which explains the various files.
   b.	We solve for a Bewley model, using the EGM method and then provide procedures to aggregate models in the space of history, and then we write a Dynare code to solve the dynamics with perturbation method.
 3.	The folder “2. BKM” provides the solution of the same model using the Boppart, Krusell, Mitman method. The main file is “MIT\_simul.jl”
 4.	The folder “3. Reiter provides the solution of the same model using the Reiter model. The main file is “Main\_Reiter.jl”. Note that we use the Dynare solver to implement the Reiter solution
 5.	The folder “4. RA” provides the solution of the same model for the representative agent economy, using Dynare.
 6.	The folder “5. Figures” contains a Julia code to Plot figures with the information contained in the other folders.

