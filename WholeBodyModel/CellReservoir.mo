model LockedPool "non-transferable cell metabolite pool"
  SingleMetabolite metabolites[MLocked] "interface for metabolites";
    
equation
  der(metabolites.M) = metabolites.Q;
  

end LockedPool;
