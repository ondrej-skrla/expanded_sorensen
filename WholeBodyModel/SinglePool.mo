model SinglePool "pool for single metabolite"
  SingleMetabolite metabolite "interface for metabolites";
  Real metabolite_mass = metabolite.M;
  
equation
  der(metabolite.M) = metabolite.Q;
  
end SinglePool;