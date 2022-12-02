model ExtraPool "extracellular metabolite pool"
  SingleMetabolite metabolites[MPlasma] "interface for metabolites";
equation
  der(metabolites.M) = metabolites.Q;
end ExtraPool;
