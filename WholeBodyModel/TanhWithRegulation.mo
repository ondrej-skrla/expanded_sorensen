model TanhWithRegulation
  extends SinkToSource;
   
  parameter Real k_1 "max reaction rate [mg/min]";
  parameter Real k_2 "steepness [-]";
  parameter Real basal_inp "basal mass of input [mg]";
  Real reg_eff "the effect of regulation [-]";
  
protected
  Real inp_normed = met_inp.M/basal_inp "normalized input [-]";
  
equation

  60*met_inp.Q = reg_eff*k_1*tanh(k_2*inp_normed);
  

end TanhWithRegulation;
