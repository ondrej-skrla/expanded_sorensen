model VariableSource
  SingleMetabolite metabolite;
  Real k_source "[mg/min]";
  
equation
  60*metabolite.Q = -k_source;

end VariableSource;