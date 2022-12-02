model TransferPool "intracellular transferable metabolite pool"
  SingleMetabolite metabolites[MChange] "interface for metabolites";
equation
  der(metabolites.M) = metabolites.Q;
end TransferPool;
