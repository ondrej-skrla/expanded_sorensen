model TanhInhibWithRegulation
  extends SinkToSource;
   
  parameter Real k_1 "reaction rate [mg/min]";
  parameter Real k_2 "steepness of inhibition [-]";
  parameter Real k_3 "half max inhibition (tanh shift) [-]";
  parameter Real basal_inp "basal mass of input [mg]";
  parameter Real basal_out "basal mass of output [mg]";
  Real reg_eff "the effect of regulation [-]";
  
protected
  Real inp_normed = met_inp.M/basal_inp "normalized input [-]";
  Real out_normed = met_out.M/basal_out "normalized output [-]";
  
equation

  60*met_inp.Q = reg_eff*k_1*inp_normed*(1+tanh(k_2*(k_3-out_normed)));
  
end TanhInhibWithRegulation;
