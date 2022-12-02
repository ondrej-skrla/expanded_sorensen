model MichealisMentenWithReg
  extends SinkToSource;
  parameter Real v_max "Maximum speed of the reaction [mg/min]";
  parameter Real k_m "Half saturation constant [-]";
  parameter Real basal_inp "basal mass of input [mg]";
  Real reg_eff "the effect of regulation [-]";
  
protected
  Real inp_normed = met_inp.M/basal_inp "normalized input"; 
  
equation
  60 * met_inp.Q = reg_eff * v_max * inp_normed / (k_m + inp_normed);

end MichealisMentenWithReg;