We provide all R source codes used for simulations

System information for testing our simulations.
1 Operation system:Linux
2 R version: 3.6.1 or newer
3 Dependencies: pcalg,CePa,reshape2,data.table,dplyr,coop,parallel,graph,grid,Rgraphviz,Corbi,glmnet,MLmetrics,MASS,igraph,stats,RNOmni,ParallelPC,glasso,glassoFast,RBGL,epiR,RiskPortfolios,matrixcalc,CVglasso,iterators,foreach,doParallel,rlist,Matrix. 
4 Main code: Simulation_Overall_Rundata.R
Note: 1)Please change the root path to your own local path to successfully load functions under folders "Simulation_Extension_Functions" and    "Simulation_Steps_Functions"
      2)Since our code includes a wide range of comparisons, it would be a bit time-consuming to get all the simulation results . 
      3)If you just want to get the main idea of our framework, please choose to run 'Simulation_ToyExample.R' under folder "ToyExample" to save time. 

5 Toy example: we have provided a toy example under foler "ToyExample" to help illustrate our framework.
  1) R code: Simulation_ToyExample.R
  2) Input:'randomDAG_sim.Rdata','gene_index.Rdata','snp_index.Rdata'
  2) Expected output: simulation_ToyExample_res1.csv, simulation_ToyExample_res2.csv
  3) Expected running time: within 30min （on our server with CPU Intel(R) Xeon(R) CPU E5-2698 v4 @ 2.20GHz）
 

